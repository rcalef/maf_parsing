/*
 * maf_parser.c
 *
 *  Created on: Aug 2, 2014
 *      Author: calef_000
 */

#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <errno.h>

#include "mafparser.h"


int get_next_offset(maf_array_parser parser){
   if(parser->curr_block >= parser->size) return -1;
   return parser->alignment_blocks[parser->curr_block++];
}

void free_sequence(seq sequence){
   if(sequence==NULL) return;
   free(sequence->src);
   free(sequence->sequence);
   free(sequence);
}

void free_alignment_block(alignment_block aln){
   if(aln==NULL) return;
   free(aln->data);
   for(int i =0; i < aln->size; ++i){
      free_sequence(aln->sequences[i]);
   }
   free(aln->sequences);
   free(aln);
}
seq copy_sequence(seq sequence){
   if(sequence==NULL) return NULL;
   seq copy = malloc(sizeof(*copy));
   copy->src=strdup(sequence->src);
   assert(copy->src != NULL);
   copy->start = sequence->start;
   copy->size = sequence->size;
   copy->strand = sequence->strand;
   copy->srcSize = sequence->srcSize;
   copy->sequence = strdup(sequence->sequence);
   assert(copy->sequence!=NULL);
   return copy;
}
seq get_sequence(char *data){
   if(data == NULL) return NULL;
   seq new_seq = malloc(sizeof(*new_seq));
   assert(new_seq!=NULL);
   new_seq->src=NULL;
   new_seq->sequence=NULL;
   char *temp = strdup(data);
   assert(temp!=NULL);
//First part of entry, is the 's', throw that away
   char *datum =strtok(temp," \t\n");
   for(int i=2;i<8;++i){
      datum = strtok(NULL," \t\n");
      if(datum == NULL){
         fprintf(stderr,"Invalid sequence: %s\n", data);
         free_sequence(new_seq);
         return NULL;
      }
   switch (i){
//Second part is species name and contig
      case 2: 
         new_seq->src = strdup(datum);
         assert(new_seq->src != NULL);
         break;
//Third part is the start of the aligned region in the source sequence
      case 3:
        errno=0;
        unsigned long start = strtol(datum,NULL,10);   
        if(errno !=0){
           fprintf(stderr, "Invalid sequence start: %s\nIn sequence: %s\n"
             ,datum,data);
           free_sequence(new_seq);
           return NULL;
        }
       new_seq->start=start;
       break;
//Fourth is aligned sequence length
      case 4:
        errno=0;
        unsigned int size = strtol(datum,NULL,10);
        if(errno !=0){
           fprintf(stderr, "Invalid sequence start: %s\nIn sequence: %s\n"
             ,datum,data);
           free_sequence(new_seq);
           return NULL; 
        }
        new_seq->size = size;
        break;
//Fifth is strand
      case 5:
        if(datum[0] != '+' && datum[0] != '-'){
           fprintf(stderr, "Invalid strand: %s\nIn sequence: %s\n"
             ,datum,data);
           free_sequence(new_seq);
           return NULL;
        }
        new_seq->strand=datum[0];
        break;
//Sixth is size of source sequence
      case 6:
        errno=0;
        unsigned long srcSize = strtol(datum,NULL,10);
        if(errno !=0){
            fprintf(stderr, "Invalid source sequence size: %s\nIn sequence: %s\n"
             ,datum,data);
            free_sequence(new_seq);
            return NULL;
        }
        new_seq->srcSize = srcSize;
        break;
//Last is the sequence itself
     case 7:
         new_seq->sequence = strdup(datum);
         assert(new_seq->sequence !=NULL);
         break;
     default:
       printf("Default case\n");
    }
  }
  free(temp);
  return new_seq;
}


alignment_block array_next_alignment(maf_array_parser parser){
   if(parser->curr_block==parser->size) return NULL;
   ++parser->curr_block;
   int check= fseek(parser->maf_file,
        parser->alignment_blocks[parser->curr_block],SEEK_SET);
   if(check !=0){
      fprintf(stderr,"File seek error: %s\n",strerror(errno));
      return NULL;
   }
   char buffer[4096];
//First read 'a' line and initialize alignment struct
   alignment_block new_align=malloc(sizeof(*new_align));
   assert(new_align != NULL);
   new_align->sequences = malloc(2*sizeof(*new_align->sequences));
   assert(new_align->sequences != NULL);
   new_align->size=new_align->curr_seq=0;
   new_align->max=2;
   char *fc = fgets(buffer,4096,parser->maf_file);
   if(ferror(parser->maf_file) != 0){
          fprintf(stderr, "File stream error: %s\nError: %s",
             parser->filename,strerror(errno));
          return NULL;
   }
//***HANDLE SCORE/PASS/DATA here***
   new_align->data = NULL;
   seq new_seq;
   while(!feof(parser->maf_file)){
      char *fc = fgets(buffer,4096,parser->maf_file);
      if(ferror(parser->maf_file) != 0){
          fprintf(stderr, "File stream error: %s\nError: %s",
             parser->filename,strerror(errno));
          return NULL;
      }
      if(buffer[0]!='s')break;
      new_seq = get_sequence(buffer);
      if(new_seq == NULL){
        fprintf(stderr, "Invalid sequence entry %s\n",buffer);
        return NULL;
      }
      if(new_align->size ==new_align->max){
          new_align->sequences=realloc(new_align->sequences,
             2*new_align->max*sizeof(seq));
          assert(new_align->sequences!=NULL);
          new_align->max *=2;
      }new_align->sequences[new_align->size++]=new_seq;
   }
   return new_align;
}

alignment_block linear_next_alignment_buffer(maf_linear_parser parser){
   alignment_block new_align = NULL;
   int in_block=0;

   long bytesread; 
   int sizeLeftover=0; 
   int bLoopCompleted = 0; 
   char *datum;
   int init=1;
   char *npos;
   do{
      if(parser->fill_buf){
         bytesread = fread(parser->buf+sizeLeftover, 1, 
            sizeof(parser->buf)-1-sizeLeftover, parser->maf_file);
         if (bytesread<1){ 
            bLoopCompleted = 1; 
            bytesread  = 0;
            continue; 
         }       
        if(ferror(parser->maf_file) != 0){
                fprintf(stderr, "File stream error: %s\nError: %s",
                   parser->filename,strerror(errno));
                return NULL;
         }
         parser->buf[sizeLeftover+bytesread]=0;
         parser->curr_pos=0;
         parser->pos=parser->buf;
         --parser->fill_buf;
         init=1;
     }
//    if(init){
//        datum= strtok(parser->buf+parser->curr_pos,"\n");
//        --init;
//     }
//     else datum = strtok(NULL,"\n");
     npos = strchr(parser->pos,'\n');
     
     if(npos==NULL){
//        datum = strtok(parser->buf+parser->curr_pos,"\n");
        sizeLeftover = strlen(parser->pos);
        memmove(parser->buf,parser->buf+(sizeof(parser->buf))-sizeLeftover-1,sizeLeftover);
        ++parser->fill_buf;
        continue;
     }
     *npos=0;
     datum = parser->pos;
     parser->pos = npos+1;
//parser->curr_pos+=strlen(datum)+1;
//If we've yet to enter an alignment block, and the first character
//of the line isn't 'a', then skip over it.
      if(!in_block && datum[0]!='a') continue;
//***HANDLE SCORE/PASS/DATA here**i
      else if(datum[0]=='a'){
//If we find an 'a' after entering a block, then this is a new block
//so rewind the file pointer and break out of read loop.
         if(in_block){
            *npos='\n';
//            *(parser->pos - 1) = '\n';
//            parser->curr_pos -= strlen(datum)-1;
            parser->pos = datum;
            break;
         }
//Else we're starting a new alignment block, initialize the data
//structure and set in_block to true.
         new_align=malloc(sizeof(*new_align));
         assert(new_align != NULL);
         new_align->sequences = malloc(16*sizeof(*new_align->sequences));
         assert(new_align->sequences != NULL);
         new_align->size=new_align->curr_seq=0;
         new_align->max=16;
         new_align->data = NULL;
         in_block=1;
         continue;
      }
//If in a block and find 's', then it's a sequence to add to the
//current alignment block, parse it, reallocate alignment block's
//sequence array if necessary, and store the new sequence.
      else if(datum[0]=='s'){
         seq new_seq = get_sequence(datum);
         if(new_seq == NULL){
           fprintf(stderr, "Invalid sequence entry %s\n",datum);
           return NULL;
         }
         if(new_align->size ==new_align->max){
             new_align->sequences=realloc(new_align->sequences,
                2*new_align->max*sizeof(seq));
             assert(new_align->sequences!=NULL);
             new_align->max *=2;
         }new_align->sequences[new_align->size++]=new_seq;
         continue;
      }
//If we hit a character other than 'a' or 's', then we've exited
//the current alignment block, break out of the read loop and return
//the current alignment block.
      else break;
   }while(!bLoopCompleted);
   return new_align;
}

alignment_block linear_next_alignment(maf_linear_parser parser){
   char buffer[4096];
   alignment_block new_align = NULL;
   int in_block=0;
   int file_pos;
   while(!feof(parser->maf_file)){
      file_pos = ftell(parser->maf_file);
      char *fc = fgets(buffer,4096,parser->maf_file);
      if(ferror(parser->maf_file) != 0){
             fprintf(stderr, "File stream error: %s\nError: %s",
                parser->filename,strerror(errno));
             return NULL;
      }
//If we've yet to enter an alignment block, and the first character
//of the line isn't 'a', then skip over it.
      if(!in_block && buffer[0]!='a') continue;
//***HANDLE SCORE/PASS/DATA here**i
      else if(buffer[0]=='a'){
//If we find an 'a' after entering a block, then this is a new block
//so rewind the file pointer and break out of read loop.
         if(in_block){
            int check= fseek(parser->maf_file,file_pos,SEEK_SET);
            if(check !=0){
               fprintf(stderr,"File seek error: %s\n",strerror(errno));
               return NULL;
            }
            break;
         }
//Else we're starting a new alignment block, initialize the data
//structure and set in_block to true.
         new_align=malloc(sizeof(*new_align));
         assert(new_align != NULL);
         new_align->sequences = malloc(16*sizeof(*new_align->sequences));
         assert(new_align->sequences != NULL);
         new_align->size=new_align->curr_seq=0;
         new_align->max=16;
         new_align->data = NULL;
         in_block=1;
         continue;
      }
//If in a block and find 's', then it's a sequence to add to the
//current alignment block, parse it, reallocate alignment block's
//sequence array if necessary, and store the new sequence.
      else if(buffer[0]=='s'){
         seq new_seq = get_sequence(buffer);
         if(new_seq == NULL){
           fprintf(stderr, "Invalid sequence entry %s\n",buffer);
           return NULL;
         }
         if(new_align->size ==new_align->max){
             new_align->sequences=realloc(new_align->sequences,
                2*new_align->max*sizeof(seq));
             assert(new_align->sequences!=NULL);
             new_align->max *=2;
         }new_align->sequences[new_align->size++]=new_seq;
         continue;
      }
//If we hit a character other than 'a' or 's', then we've exited
//the current alignment block, break out of the read loop and return
//the current alignment block.
      else break;
   }
   return new_align;
}



void array_double(maf_array_parser parser){
   parser->alignment_blocks =realloc(parser->alignment_blocks,2*parser->max*sizeof(int));
   assert(parser->alignment_blocks!=NULL);
   parser->max *=2;
   printf("Doubling array to max size of %d\n", parser->max);
}
maf_linear_parser get_linear_parser(FILE *maf_file, char *filename){
	maf_linear_parser parser = malloc(sizeof(*parser));
	assert(parser!=NULL);
	parser->maf_file = maf_file;
	parser->filename= strdup(filename);
	assert(filename!=NULL);
        parser->curr_pos=0;
        parser->fill_buf=1;
	return parser;
}

maf_array_parser get_array_parser(FILE *maf_file,char *filename){
	char buffer[3000];
	assert(buffer != NULL);
	maf_array_parser parser = malloc(sizeof(*parser));
	assert(parser != NULL);
        parser->maf_file = maf_file;
        parser->filename = strdup(filename);
        assert(parser->filename != NULL);
        parser->curr_block=0;
        parser->size = 0;
        parser->alignment_blocks = malloc(2*sizeof(int));
        assert(parser->alignment_blocks != NULL);
        parser->max = 2;
	int pos;
	while(!feof(maf_file)){
		pos = ftell(maf_file);
		char *check = fgets(buffer,3000,maf_file);
		if(ferror(maf_file) != 0){
			fprintf(stderr, "File stream error: %s\nError: %s",
					filename,strerror(errno));
			return NULL;
		}
		if(buffer[0]=='a'){
			if(parser->size == parser->max)
                           array_double(parser);
                        parser->alignment_blocks[parser->size++]=pos;
		}
	}
	return parser;
}

seq iterate_sequences(alignment_block aln){
   if(++aln->curr_seq ==aln->size) return NULL;
   return aln->sequences[aln->curr_seq];
}

void free_linear_parser(maf_linear_parser parser){
   free(parser->filename);
   free(parser);
}

void free_array_parser(maf_array_parser parser){
    free(parser->filename);
    free(parser->alignment_blocks);
    free(parser);
    return;
}
void print_sequence(seq sequence){
   if(sequence==NULL) return;
   printf("s %s  %lu  %u  %c  %lu  %s\n"
      ,sequence->src,sequence->start,sequence->size
      ,sequence->strand,sequence->srcSize,sequence->sequence);
}
void print_alignment(alignment_block aln){
   if(aln==NULL)return;
   printf("\na %s\n",aln->data);
   for(int i=0; i < aln->size; ++i)
    if(aln->sequences[i] != NULL) print_sequence(aln->sequences[i]);
}