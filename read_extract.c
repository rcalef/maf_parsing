#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>

//Note that hash table functions and data types are included via mafparser.h,
//which includes the standard C library 'search.h'
#include "mafparser.h"

typedef struct _read_positions{
   char *scaf_name;
   char *read_name;
   char *aligned_seq;
   int *positions;
   unsigned int num_positions;
   unsigned int max_positions;
}*read_positions;

//Next define global variables used to keep track of parameters passed in to
char **tree_species_names;
int tree_species_size;
int tree_species_max;

//Define long options, note that options with 'no_argument'
//specified do require arguments, the 'no_argument' specification
//is necessary for reading in a variable size list of arguments.
static struct option long_options[]={
{"tree-species",no_argument,0,'t'},
//Ignore the line below, necessary to define the long options
{0,0,0,0}
};

void free_read(read_positions read){
   free(read->scaf_name);
   free(read->read_name);
   free(read->positions);
   free(read->aligned_seq);
   free(read);
}
/*
get_largest_index takes in a list of numbers, and the size of the list
and returns the index of the largest element in the size, or -1 if the
list is empty (somewhat poor design choice, but in this program,
this function will only take in lists of integers greater than or 
equal to zero)
*/
int get_largest_index(int *nums, int size){
   int max = 0;
   int index = -1;
   for(int i =0; i < size; ++i){ 
      if(nums[i] > max){
         max = nums[i];
         index=i;
      }
   }
   return index;
}

/*
search_hash is a convenience function for checking whether or not the given
hash table contains an entry for the key of interest.

Takes in the hash table to search, the key to search for, and a pointer to
an empty ENTRY object. The hsearch method will either store the found entry
in 'ret_val', or return 'ret_val' as NULL if the key is not found.
*/
ENTRY *search_hash(char *key, ENTRY *ret_val, hash table){
         ENTRY search={key,NULL};
         int hc = hsearch_r(search,FIND,&ret_val,table);
         if(hc == 0){
            //errno == 3 corresponds to key not found, hence return NULL,
            //otherwise some other, real error occurred.
            if(errno == 3) return NULL;
            fprintf(stderr,"Error searching hash table: %s\n", strerror(errno));
            exit(1);
         }
         return ret_val;
}

//Convenience function to free up memory when the program exits.
void clean_up(){
   free(tree_species_names);
}

/*
parse_args is the convenience function gather together all command line
argument parsing, storing the arguments in the appropriate global variables
defined at the beginning of the file
*/
void parse_args(int argc, char **argv){
   //Minimum of 3 arguments, in group, out group, and MAF file
   if (argc < 2){
      fprintf(stderr,"Too few arguments\n");
      exit(1);
   }
   char c;
   int option_index=0;
   while((c=getopt_long(argc,argv,"t",long_options,&option_index))!= -1){
      switch(c){
         //'-g' option to specify which species to output an annotated genome
         //for, parsed in the same manner as -i and -o
         case 't':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "--tree-species parameter requires at least one argument\n");
               exit(1);
            }
            do{
               if(strcasestr(argv[optind],".bed")!=NULL) return;
               if(tree_species_size == tree_species_max){
                  tree_species_max *=2;
                  tree_species_names=realloc(tree_species_names,
                     tree_species_max*sizeof(char*));
               }
               tree_species_names[tree_species_size++]=argv[optind++];
            }while(optind < argc && argv[optind][0]!='-');
            break;
         //Default case to catch invalid options
         case '?':
	    if(optopt == NULL) fprintf(stderr,"Invalid long option: %s\n",argv[optind-1]);
	    else fprintf(stderr, "Invalid short option: %c\n", optopt);
            exit(1);
      }
   }
}


read_positions cigar_to_positions(char *cigar, unsigned int start,char *sequence){
   unsigned int i = 0;
   unsigned int num_bases = 0;
   unsigned int sequence_pos = 0;
   unsigned int num_positions = 0;
   unsigned int max_positions = 128;
   int *positions = malloc(sizeof(int) * max_positions);
   assert(positions != NULL);
   char *align_seq = calloc(max_positions,sizeof(char));
   assert(align_seq != NULL);
   char *tmp = NULL;
   char tmp_char;
   while(1){
      tmp = strpbrk(cigar,"MIDN");
      if(tmp == NULL) break;
      tmp_char = *tmp;
      *tmp = '\0';
      num_bases = atoi(cigar);
      switch(tmp_char){
         case 'M':
            for(i = start; i < start + num_bases; ++i){
               if(num_positions >= max_positions){
                  max_positions *=2;
                  positions = realloc(positions,
                                sizeof(unsigned int) * max_positions);
                  assert(positions != NULL);
                  align_seq = realloc(align_seq,
                                sizeof(char) * max_positions);
                  assert(align_seq != NULL);
               }
               align_seq[num_positions] = sequence[sequence_pos++];
               positions[num_positions++] = i;
            }
            start = i;
            break;
         case 'I':
            for(i = 0; i < num_bases; ++i){
               if(num_positions >= max_positions){
                  max_positions *=2;
                  positions = realloc(positions,
                                sizeof(unsigned int) * max_positions);
                  assert(positions != NULL);
                  align_seq = realloc(align_seq,
                                sizeof(char) * max_positions);
                  assert(align_seq != NULL);

               }
               align_seq[num_positions] = sequence[sequence_pos++];
               positions[num_positions++] = -1;
            }
            break;
         case 'D':
            for(i = 0; i < num_bases; ++i){
               if(num_positions >= max_positions){
                  max_positions *=2;
                  positions = realloc(positions,
                                sizeof(unsigned int) * max_positions);
                  assert(positions != NULL);
                  align_seq = realloc(align_seq,
                                sizeof(char) * max_positions);
                  assert(align_seq != NULL);

               }
               align_seq[num_positions] = '-';
               positions[num_positions++] = start + i;
            }
            start += num_bases;
            break;
         default:
            fprintf(stderr,"Unexpected character in cigar string: %c\n",
               tmp_char);
            break;
      }
      if(*(tmp + 1) == '\0') break;
      cigar = tmp + 1;
   }
   read_positions new_pos = malloc(sizeof(struct _read_positions));
   assert(new_pos != NULL);
   align_seq[num_positions] = '\0';
   new_pos->positions = positions;
   new_pos->num_positions = num_positions;
   new_pos->max_positions = max_positions;
   new_pos->aligned_seq = align_seq;
   return new_pos;
}




void single_read_region(read_positions read){
   int i = 0;
   char filename[128];
   maf_linear_parser scaf_parser = NULL;
   sorted_alignment_block aln = NULL;
   seq reference_seq = NULL;
   unsigned int current_column = 0;
   unsigned int aln_position = 0;
   unsigned int ref_position = 0;
   unsigned int pos_index = 0;
   unsigned int end_range = 0;
   char **multiple_aln = malloc(sizeof(char *) * tree_species_size);
   assert(multiple_aln != NULL);
   for(i = 0; i < tree_species_size; ++i){
      multiple_aln[i] = malloc(sizeof(char) * (read->num_positions+1));
      assert(multiple_aln[i] != NULL);
      memset(multiple_aln[i],45,read->num_positions);
      memset(multiple_aln[i] + read->num_positions,0,1);
   }
   FILE *curr_scaf_maf = NULL;
   FILE *fasta_output = NULL;
   strcpy(filename,read->read_name);
   strcat(filename,"-regions.fasta");
   if((fasta_output= fopen(filename, "w")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s\n",
         filename,strerror(errno));
      exit(1);
   }
   strcpy(filename,"/projects/redser/ftp/100way_alignment/split_maf/birdRepAnc05.");
   strcat(filename,read->scaf_name);
   strcat(filename,"-blocks_only.maf");
   if((curr_scaf_maf = fopen(filename,"r")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s\n",
         filename,strerror(errno));
      exit(1);
   }
   scaf_parser = get_linear_parser(curr_scaf_maf,filename);
   while(1){
      aln = get_sorted_alignment(scaf_parser,tree_species_names,
         tree_species_size, NULL,0);
      if(aln == NULL) break;
      aln_position = 0;
      for(int j = 0; j < aln->in_size; ++j){
//           if(strcmp(aln->in_sequences[j]->species,"birdRepAnc05") == 0){
         if(strcmp(aln->in_sequences[j]->species,"birdRepAnc05") == 0){
            reference_seq = aln->in_sequences[j];
            break;
         }
      }
      ref_position = reference_seq->start;
      end_range = ref_position + reference_seq->size;
      if(read->positions[pos_index] >= end_range){
         free_sorted_alignment(aln);
         continue;
      }
      print_sorted_alignment(aln);
      while(1){
         if(pos_index >= read->num_positions){
            free_sorted_alignment(aln);
            goto done;
         }
         if(read->positions[pos_index] < ref_position){
            ++pos_index;
            ++current_column;
            continue;
         }
         //-1 corresponds to I, insertion in read relative to reference, so
         //need to add gaps to MAF sequences. As all sequences were initialized
         //to a string of gaps, we can just increment the column and skip.
         if(read->positions[pos_index] == -1){
            ++pos_index;
            ++current_column;   
            continue;
         }
         //-2 corresponds to D, deletion in read relative to reference, gaps
         //are already present in read->aligned seq, so just add the correct
         //bases to the species sequences as normal.
//         if(positions[pos_index] == -2){
//            ++pos_index;
//            ++cirr
//
//         }
         if(read->positions[pos_index] >= end_range) break;
         for(int j = aln_position; j < aln->seq_length; ++j){
            if(reference_seq->sequence[j] != '-'){
               if(ref_position == read->positions[pos_index]){
                  aln_position = j;
                  break;
               }
               ++ref_position;
            }
         }
         fprintf(stderr,"Position %u in scaffold %s found in alignment column %d\n",
            read->positions[pos_index],read->scaf_name,aln_position);
         for(int j = 0; j < aln->in_size; ++j){
            for(int k = 0; k < tree_species_size; ++k){
               if(strcmp(aln->in_sequences[j]->species,
                  tree_species_names[k]) !=0) continue;
//                  if(strcmp(aln->in_sequences[j]->scaffold,
//                     curr_scaf->scaf_name) != 0) continue;
               multiple_aln[k][current_column] = aln->in_sequences[j]->sequence[aln_position];
               fprintf(stderr,"%s\t%c\n",tree_species_names[k],aln->in_sequences[j]->sequence[aln_position]);
            }
         }
         ++pos_index;
         ++current_column;
      }
      free_sorted_alignment(aln);
      if(pos_index >= read->num_positions) break;
   }
   done:
   fprintf(fasta_output,">%s\n%s\n",read->read_name,read->aligned_seq);
   for(i = 0; i < tree_species_size; ++i){
      fprintf(fasta_output,">%s\n%s\n", tree_species_names[i],multiple_aln[i]);
      free(multiple_aln[i]);
   }
   free(multiple_aln);
   free_linear_parser(scaf_parser);
   fclose(curr_scaf_maf);
   fclose(fasta_output);
}

void process_bed_complete(FILE *bed_file){
   char buffer[1024];
   char *scaf_name = NULL;
   char *read_name = NULL;
   char *cigar = NULL;
   char *sequence = NULL;
   read_positions curr_read = NULL;
   unsigned int i = 0;
   unsigned int start = 0;
   unsigned int end = 0;
   while(fgets(buffer,1024,bed_file) != NULL){
      scaf_name = strtok(buffer,"\t");
      start = atoi(strtok(NULL,"\t"));
      end = atoi(strtok(NULL,"\t"));
      read_name = strtok(NULL,"\t");
      sequence = strtok(NULL,"\t");
      cigar = strtok(NULL,"\t");
      cigar = strtok(NULL,"\t");
      curr_read = cigar_to_positions(cigar,start,sequence);
      curr_read->read_name = strdup(read_name);
      assert(curr_read->read_name != NULL);
      curr_read->scaf_name = strdup(scaf_name);
      assert(curr_read->scaf_name != NULL);
      single_read_region(curr_read);
      free_read(curr_read);
   }
}



void set_up(){
   tree_species_size=0;
   tree_species_max=8;
   tree_species_names = malloc(sizeof(char *)*tree_species_max);
}

void print_multiple_aln(char ** multiple_aln){
   for(int i = 0; i < tree_species_size; ++i)
      printf(">%s\n%s\n", tree_species_names[i],multiple_aln[i]);
}

int main(int argc, char **argv){
   set_up();
   parse_args(argc,argv);
   if(optind >= argc){
      fprintf(stderr, "Missing required BED filename\n");
      exit(1);
   }
   char *filename = argv[optind];
   FILE *bed_file;
   if((bed_file= fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s\n",
         filename,strerror(errno));
      return 1;
   }
   int i = 0;
   for(i = 0; i < tree_species_size; ++i) fprintf(stderr,"%s\n", tree_species_names[i]);
   fprintf(stderr,"Filename: %s\n",filename);
   fprintf(stderr,"Num tree species: %d\n",tree_species_size);
   process_bed_complete(bed_file);
   fclose(bed_file);
   clean_up();
   return 0;
}
