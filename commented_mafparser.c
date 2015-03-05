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

/*
in_list simply returns true if the string 'needle' is in the unordered
array of strings 'haystack' which has 'size' elements.
*/
int in_list(char *needle, char **haystack, int size){
   for(int i = 0; i < size; ++i)
      if(!strcmp(needle,haystack[i])) return 1;
   return 0;
}



/*
free_sequence cleans up the memory used by a sequence struct
*/
void free_sequence(seq sequence){
   if(sequence==NULL) return;
   free(sequence->src);
   if(sequence->sequence != NULL)
     free(sequence->sequence);
   free(sequence->species);
   free(sequence);
}

/*
free_sorted_alignment cleans up the memory used by a sorted_alignment_block 
struct
*/
void free_sorted_alignment(sorted_alignment_block aln){
   if(aln==NULL) return;
   free(aln->data);
   int i =0;
   for(; i < aln->in_size; ++i){
      free_sequence(aln->in_sequences[i]);
   }
   free(aln->in_sequences);
   for(i=0; i < aln->out_size; ++i){
      free_sequence(aln->out_sequences[i]);
   }
   free(aln->out_sequences);
   free(aln);
}

/*
Convenience function to make a new sequence struct as a copy of an existing
one, i.e. a copy constructor
*/
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
   copy->species = strdup(sequence->species);
   assert(copy->species != NULL);
   copy->scaffold = strdup(sequence->scaffold);
   assert(copy->scaffold != NULL);
   return copy;
}

/*
Takes in a string containing a MAF alignment block sequence entry and parses
the data in to a new sequence struct.
*/
seq get_sequence(char *data){
   if(data == NULL) return NULL;
   char *seq_parse;
   char *src_parse;
   seq new_seq = malloc(sizeof(*new_seq));
   assert(new_seq!=NULL);
   new_seq->src=NULL;
   new_seq->sequence=NULL;
   char *temp = strdup(data);
   assert(temp!=NULL);
   //First part of entry, is the 's', throw that away
   char *datum =strtok_r(temp," \t\n",&seq_parse);
   for(int i=2;i<8;++i){
      datum = strtok_r(NULL," \t\n",&seq_parse);
      if(datum == NULL){
         fprintf(stderr,"Invalid sequence: %s\n", data);
         free_sequence(new_seq);
         return NULL;
      }
      switch (i){
         //Second part is species name and contig
         case 2: 
            new_seq->src = strdup(datum);
            char *parse_src = strdup(datum);
            new_seq->species = strtok_r(parse_src,".",&src_parse);
            new_seq->scaffold = strtok_r(NULL,".",&src_parse);
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



/*
sorted_alignment_block takes in an existing maf_linaer_parser struct, lists of
in-group and out-group species, and their respective sizes. 

The function returns the next alignment block in the MAF file, with sequence 
entries split in to two lists, one containing entries for in-group species and 
one with entries for out-group species. Sequence entries for species in neither
group will be thrown away.
*/
sorted_alignment_block get_sorted_alignment(maf_linear_parser parser, 
              char **in_group, int in_size, char **out_group, int out_size){
   //First initialize all the variables used for the parsing loop: the 
   //allignment block object, 'in_block' to indicate whether we're reading
   //lines of the alignment block.
   //
   //The core of the parsing loop uses a large buffer, reading data from the
   //file in large chunks of bytes. Use 'bytesread' and'sizeLeftover' to keep
   //track of the number of bytes being read in and kept between buffer
   //fillings.
   //
   //In general, the parser contains a buffer that store a large chunk of data
   //with no information as to where it ends, so it probably ends in the
   //middle of a line. To handle this problem, we do multiple things:
   //   First, always keep track of the last new line character seen 
   //   Next, at the end of the buffer, move the bytes since the last newline
   //    to the beginning of the buffer, overwriting whatever old data:
   //       memmove(parser->buf,parser->buf+(sizeof(parser->buf))-sizeLeftover-1
   //             ,sizeLeftover);
   //   Third, read in new data starting right after the leftover data ends:
   //       bytesread = fread(parser->buf+sizeLeftover, 1,
   //             sizeof(parser->buf)-1-sizeLeftover, parser->maf_file);
   //     and then we can continue parsing line by line.
   sorted_alignment_block new_align = NULL;
   int in_block=0;
   long bytesread;
   int sizeLeftover=0;
   int bLoopCompleted = 0;
   int first = 1;
   char *datum;
   char *npos;
   char *temp;
   char *species;
   do{
      //The parser's fill_buf variable indicates whether or not the parser's
      //buffer needs to be refilled.
      if(parser->fill_buf){
         /Store the number of bytes read, primarily to check if we've read the
         //whole file. The minus one in the fread call allows room for the null
         //plug at the end of the buffer, making the buffer a legal C string.
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
         //Write a null plug to the end of the buffer, and reset positions
         parser->buf[sizeLeftover+bytesread]=0;
         parser->curr_pos=0;
         parser->pos=parser->buf;
         --parser->fill_buf;
      }

      //In the general case of the loop, we just read a line, and parse it.
      //First get the position of the next new line character, if none are found
      //then we need to refill the buffer.
      npos = strchr(parser->pos,'\n');
      if(npos==NULL){
         sizeLeftover = strlen(parser->pos);
         memmove(parser->buf,parser->buf+(sizeof(parser->buf))-sizeLeftover-1,
            sizeLeftover);
         ++parser->fill_buf;
         continue;
      }
    
      //If the whole line is there, change the new line character to a null plug
      //so the individual line is now a C string, and adjust the parser position
      *npos=0;
      datum = parser->pos;
      parser->pos = npos+1;

      //If we've yet to enter an alignment block, and the first character
      //of the line isn't 'a', then skip over it ('a' line indicates the start
      //of a new alignment block).
      if(!in_block && datum[0]!='a') continue;
//***HANDLE SCORE/PASS/DATA here**
      //If we find a line beginning with 'a' after entering a block, then 
      //this is a new block and we're done parsiong, so rewind the parser
      //position and reset the new line character.
      else if(datum[0]=='a'){
         //If we find a line beginning with 'a' after entering a block, then 
         //this is a new block and we're done parsing, so rewind the parser
         //position and reset the new line character.
         if(in_block){
            *npos='\n';
            parser->pos = datum;
            break;
         }
         //Else we're starting a new alignment block, initialize the data
         //structure and set in_block to true.
         new_align=malloc(sizeof(*new_align));
         assert(new_align != NULL);
         new_align->in_sequences = malloc(16*sizeof(*new_align->in_sequences));
         assert(new_align->in_sequences != NULL);
         new_align->in_size=0;
         new_align->in_max=16;
         new_align->out_sequences = malloc(16*sizeof(*new_align->out_sequences));
         assert(new_align->out_sequences != NULL);
         new_align->out_size=0;
         new_align->out_max=16;
         new_align->data = NULL;
         new_align->seq_length=0;
         in_block=1;
         continue;
      }
      //If in a block and find 's', then it's a sequence to add to the
      //current alignment block, parse it, check if its species is in the
      //in- our out-groups, if so, reallocate alignment block's
      //sequence array if necessary, and store the new sequence.
      else if(datum[0]=='s'){
         seq new_seq = get_sequence(datum);
         if(new_seq == NULL){
            fprintf(stderr, "Invalid sequence entry %s\n",datum);
            return NULL;
         }
         //As of right now, making the kludgy fix of using the 
//********************************************************************
         if(first){
            new_align->seq_length=strlen(new_seq->sequence);
            first = 0;
         }
         temp = strdup(new_seq->src);
         assert(temp != NULL);
         species = strtok(temp,".");
         if(in_list(species,in_group,in_size)){
            if(new_align->in_size ==new_align->in_max){
               new_align->in_sequences=realloc(new_align->in_sequences,
                  2*new_align->in_max*sizeof(seq));
               assert(new_align->in_sequences!=NULL);
               new_align->in_max *=2;
            }new_align->in_sequences[new_align->in_size++]=new_seq;
         }else if(in_list(species,out_group,out_size)){
            if(new_align->out_size ==new_align->out_max){
               new_align->out_sequences=realloc(new_align->out_sequences,
                  2*new_align->out_max*sizeof(seq));
               assert(new_align->out_sequences!=NULL);
               new_align->out_max *=2;
            }new_align->out_sequences[new_align->out_size++]=new_seq;
         //If not in in group or out group, throw away.
         }else free_sequence(new_seq);
         free(temp);
      }
      //If we hit a character other than 'a' or 's', then we've exited
      //the current alignment block, break out of the read loop and return
      //the current alignment block.
      else break;
   }while(!bLoopCompleted);
   return new_align;
}

/*
get_linear_parser is a constructor for a new maf_linear_parser object,
taking in a FILE object and filename.
*/
maf_linear_parser get_linear_parser(FILE *maf_file, char *filename){
	maf_linear_parser parser = malloc(sizeof(*parser));
	assert(parser!=NULL);
	parser->maf_file = maf_file;
	parser->filename= strdup(filename);
	assert(filename!=NULL);
        parser->curr_pos=0;
        parser->fill_buf=1;
	parser->pos=0;
	return parser;
}

/*
Destructor for a maf_linear_parser struct
DOES NOT CLOSE THE FILE
*/
void free_linear_parser(maf_linear_parser parser){
   free(parser->filename);
   free(parser);
}

/*
Convenience function for printing a sequence object to stdout
*/
void print_sequence(seq sequence){
   if(sequence==NULL) return;
   printf("s %25s  %18lu  %8u  %c  %18lu  %s\n"
      ,sequence->src,sequence->start,sequence->size
      ,sequence->strand,sequence->srcSize,sequence->sequence);
}

/*
Convenience function for printing a sorted_alignment_block object to stdout
*/
void print_sorted_alignment(sorted_alignment_block aln){
   if(aln==NULL) return;
   printf("\na %s\n", aln->data);
   int i =0;
   for(;i<aln->in_size; ++i) print_sequence(aln->in_sequences[i]);
   for(i=0; i < aln->out_size; ++i) print_sequence(aln->out_sequences[i]);
}

