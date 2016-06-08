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

typedef struct _scaffold_positions{
   char *scaf_name;
   unsigned int *positions;
   unsigned int num_positions;
   unsigned int max_positions;
}*scaf_positions;

//Next define global variables used to keep track of parameters passed in to
char **tree_species_names;
int tree_species_size;
int tree_species_max;
int num_columns;
scaf_positions *ones_positions;
int num_scafs;
int max_scafs;
char mark_char;

//Define long options, note that options with 'no_argument'
//specified do require arguments, the 'no_argument' specification
//is necessary for reading in a variable size list of arguments.
static struct option long_options[]={
{"tree-species",no_argument,0,'t'},
{"num-columns",required_argument,0,'n'},
{"character",required_argument,0,'c'},
//Ignore the line below, necessary to define the long options
{0,0,0,0}
};


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
   for(int i = 0; i < num_scafs; ++i){
      free(ones_positions[i]->scaf_name);
      free(ones_positions[i]->positions);
      free(ones_positions[i]);
   }
   free(ones_positions);
}

/*
parse_args is the convenience function gather together all command line
argument parsing, storing the arguments in the appropriate global variables
defined at the beginning of the file
*/
void parse_args(int argc, char **argv){
   //Minimum of 3 arguments, in group, out group, and MAF file
   if (argc < 4){
      fprintf(stderr,"Too few arguments\n");
      exit(1);
   }
   char c;
   int option_index=0;
   while((c=getopt_long(argc,argv,"n:c:t",long_options,&option_index))!= -1){
      switch(c){
         //'-g' option to specify which species to output an annotated genome
         //for, parsed in the same manner as -i and -o
         case 't':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "--tree-species parameter requires at least one argument\n");
               exit(1);
            }
            do{
               if(strcasestr(argv[optind],".fasta")!=NULL) return;
               if(tree_species_size == tree_species_max){
                  tree_species_max *=2;
                  tree_species_names=realloc(tree_species_names,
                     tree_species_max*sizeof(char*));
               }
               tree_species_names[tree_species_size++]=argv[optind++];
            }while(optind < argc && argv[optind][0]!='-');
            break;
         case 'n':
            if(optarg==NULL || optarg[0]=='-'){
               fprintf(stderr, "--num-columns parameter requires one argument\n");
               exit(1);
            }
            num_columns=atoi(optarg);
            if(num_columns == 0){
               fprintf(stderr, "Invalid number of alignment columns: %d\n",num_columns);
               exit(1);
            }
            break;
         case 'c':
            if(optarg==NULL || optarg[0]=='-'){
               fprintf(stderr, "--num-columns parameter requires one argument\n");
               exit(1);
            }
            if(optarg[0] != '0' && optarg[0] != '1' && optarg[0] != '2'){
               fprintf(stderr,"--character parameter must be 0, 1, or 2\n");
               exit(1);
            }
            mark_char = optarg[0];
            break;
         //Default case to catch invalid options
         case '?':
	    if(optopt == NULL) fprintf(stderr,"Invalid long option: %s\n",argv[optind-1]);
	    else fprintf(stderr, "Invalid short option: %c\n", optopt);
            exit(1);
      }
   }
}

void process_fasta(FILE *fasta_input){
   int i = 0;
   int num_ones = 0;
   char buffer[1024];
   char *scaffold_name = NULL;
   scaf_positions curr_scaf = NULL;
   unsigned int num_positions = 0;
   unsigned int max_positions = 128;
   unsigned int *positions = malloc(sizeof(unsigned int) * 128);
   assert(positions != NULL);
   unsigned int scaffold_position = 0;
   while(fgets(buffer,1024,fasta_input) != NULL){
      if(buffer[0] == '>'){
         if(curr_scaf != NULL){
            curr_scaf->num_positions = num_positions;
            curr_scaf->max_positions = max_positions;
            curr_scaf->positions = positions;
            curr_scaf->scaf_name = scaffold_name;
            num_positions = 0;
            max_positions = 128;
            positions = malloc(sizeof(unsigned int) * 128);
            assert(positions != NULL);
            if(num_scafs >= max_scafs){
               max_scafs *=2;
               ones_positions = realloc(ones_positions,
                  sizeof(struct _scaffold_positions) * max_scafs);
               assert(ones_positions != NULL);
            }
            ones_positions[num_scafs++] = curr_scaf;
            scaffold_position = 0;
         }
         scaffold_name = strtok(buffer,".\n ");
         scaffold_name = strdup(strtok(NULL,"\n "));
         assert(scaffold_name != NULL);
         curr_scaf = malloc(sizeof(struct _scaffold_positions));
         assert(curr_scaf != NULL);
         continue;
      }
      for(i = 0; i < 1024; ++i){
         if(buffer[i] == '\n') break;
         if(buffer[i] == mark_char){
            if(num_positions >= max_positions){
               max_positions *=2;
               positions = realloc(positions,
                  sizeof(unsigned int) * max_positions);
               assert(ones_positions != NULL);
            }
            positions[num_positions++] = scaffold_position + i;
            ++num_ones;
            if(num_ones >= num_columns) goto break_loop;
         } 
      }
      scaffold_position += i;       
   }
   break_loop: 
   curr_scaf->num_positions = num_positions;
   curr_scaf->max_positions = max_positions;
   curr_scaf->positions = positions;
   curr_scaf->scaf_name = scaffold_name;
   if(num_scafs >= max_scafs){
      max_scafs *=2;
      ones_positions = realloc(ones_positions,
         sizeof(struct _scaffold_positions) * max_scafs);
      assert(ones_positions != NULL);
   }
   ones_positions[num_scafs++] = curr_scaf;
}


char **get_tree_sequences(){
   int i = 0;
   char filename[128];
   FILE *curr_scaf_maf = NULL;
   maf_linear_parser scaf_parser = NULL;
   sorted_alignment_block aln = NULL;
   scaf_positions curr_scaf = NULL;
   seq reference_seq = NULL;
   int current_column = 0;
   int aln_position = 0;
   int ref_position = 0;
   unsigned int current_pos = 0;
   unsigned int end_range = 0;
   char **multiple_aln = malloc(sizeof(char *) * tree_species_size);
   assert(multiple_aln != NULL);
   for(i = 0; i < tree_species_size; ++i){
      multiple_aln[i] = malloc(sizeof(char) * (num_columns+1));
      memset(multiple_aln[i],45,num_columns);
      memset(multiple_aln[i] + num_columns,0,1);
      assert(multiple_aln[i] != NULL);
   }
   for(i = 0; i < num_scafs; ++i){
      curr_scaf = ones_positions[i];
      current_pos = 0;
//      strcpy(filename,"birdRepAnc05.");
      strcpy(filename,"/projects/redser/ftp/100way_alignment/split_maf/birdRepAnc05.");
      strcat(filename,curr_scaf->scaf_name);
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
//            if(strcmp(aln->in_sequences[j]->species,"birdRepAnc05") == 0){
            if(strcmp(aln->in_sequences[j]->species,"birdRepAnc05") == 0){
               reference_seq = aln->in_sequences[j];
               break;
            }
         }
         ref_position = reference_seq->start;
         end_range = reference_seq->start + reference_seq->size;
         while(1){
            if(current_pos >= curr_scaf->num_positions){
               free_sorted_alignment(aln);
               goto next_scaf;   
            }
            if(curr_scaf->positions[current_pos] >= end_range) break;
            for(int j = aln_position; j < aln->seq_length; ++j){
               if(reference_seq->sequence[j] != '-'){
                  if(ref_position == curr_scaf->positions[current_pos]){
                     aln_position = j;
                     break;
                  }
                  ++ref_position;
               }
            }
            fprintf(stderr,"Position %u in scaffold %s found in alignment column %d\n",
               curr_scaf->positions[current_pos],curr_scaf->scaf_name,aln_position);
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
            ++current_pos;
            ++current_column;
            if(current_column >= num_columns) break;
         }
         free_sorted_alignment(aln);
         if(current_column >= num_columns) break;
         if(current_pos >= curr_scaf->num_positions) break;
      }
      next_scaf:
      free_linear_parser(scaf_parser);
      fclose(curr_scaf_maf);
   }
   return multiple_aln; 
}


void set_up(){
   tree_species_size=0;
   tree_species_max=8;
   tree_species_names = malloc(sizeof(char *)*tree_species_max);
   num_columns = 0;
   num_scafs = 0;
   max_scafs = 32;
   ones_positions = malloc(sizeof(struct _scaffold_positions) * max_scafs);
   assert(ones_positions != NULL);
   mark_char = '1';
}

void print_multiple_aln(char ** multiple_aln){
//   for(int i = 0; i < tree_species_size; ++i)
//      printf("%s\t%s\n", tree_species_names[i],multiple_aln[i]);
   for(int i = 0; i < tree_species_size; ++i)
      printf(">%s\n%s\n", tree_species_names[i],multiple_aln[i]);
}

int main(int argc, char **argv){
   set_up();
   parse_args(argc,argv);
   if(optind >= argc){
      fprintf(stderr, "Missing required FASTA filename\n");
      exit(1);
   }
   char *filename = argv[optind];
   FILE *fasta_file;
   if((fasta_file= fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s\n",
         filename,strerror(errno));
      return 1;
   }
   int i = 0;
   for(i = 0; i < tree_species_size; ++i) fprintf(stderr,"%s\n", tree_species_names[i]);
   fprintf(stderr,"Number of alignment columns: %d\n", num_columns);
   fprintf(stderr,"Filename: %s\n",filename);
   fprintf(stderr,"Num tree species: %d\n",tree_species_size);
   process_fasta(fasta_file);
//   for(i = 0; i < num_scafs; ++i){
//      printf("%s\n",ones_positions[i]->scaf_name);
//      for(int j = 0; j < ones_positions[i]->num_positions; ++j)
//         printf("%u ",ones_positions[i]->positions[j]);
//   }
   char ** multiple_aln = get_tree_sequences();
   for(i = 0; i < tree_species_size; ++i){
      printf(">%s\n%s\n", tree_species_names[i],multiple_aln[i]);
      free(multiple_aln[i]);
   }
   free(multiple_aln);
   fclose(fasta_file);
   clean_up();
   return 0;
}
