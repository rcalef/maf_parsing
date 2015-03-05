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


//Next define global variables used to keep track of parameters passed in to
//the program. 'in_group' and 'out_group' store the species names of the in
//group and the out group, '*_size' and '*_max' store the current size and
//maximum size of each list respectively.
char **in_group;
int in_size;
int in_max;
char **out_group;
int out_size;
int out_max;
//Also store the conservation thresholds, one for conservation in the in group
//and one for conservation in the out group.
double in_cons_thresh;
double out_cons_thresh;
//Finally define a hash table of genomes, using species names as keys, along
//with variables to keep track of the current size and max size. Like above,
//keep a list of the genome names, i.e. the keys of the hash table, allowing
//iteration over the elements of the table.
char **tree_species_names;
int tree_species_size;
int tree_species_max;
int num_columns;

//Define long options, note that options with 'no_argument'
//specified do require arguments, the 'no_argument' specification
//is necessary for reading in a variable size list of arguments.
static struct option long_options[]={
{"in-thresh",required_argument,0,'x'},
{"out-thresh",required_argument,0,'z'},
{"out-group",no_argument,0,'o'},
{"in-group",no_argument,0,'i'},
{"tree-species",no_argument,0,'t'},
{"num-columns",required_argument,0,'n'},
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
   free(in_group);
   free(out_group);
//   for(int i = 0; i < tree_species_size; ++i){
//      free(tree_species_names[i]);
//   }
   free(tree_species_names);
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
   while((c=getopt_long(argc,argv,"n:x:z:iot",long_options,&option_index))!= -1){
      switch(c){
         //'-i' option to specify in group species, takes in an arbitrary sized
         //list of species
         case 'i':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "--in-group parameter requires at least one argument\n");
               exit(1);
            }
            do{
               //Need to break out of arg parsing if we hit the MAF file argument
               if(strcasestr(argv[optind],".maf")!=NULL) return;
               //Reallocate in group species list if full
               if(in_size == in_max){
                  in_max*=2;
                  in_group=realloc(in_group,
                     in_max*sizeof(char*));
               }
               in_group[in_size++]=argv[optind++];
            }while(optind < argc && argv[optind][0]!='-');
            break;
         //'-o' option to specify out group species, just like above for in group
         case 'o':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "--out-group parameter requires at least one argument\n");
               exit(1);
            }
            do{
               if(strcasestr(argv[optind],".maf")!=NULL) return;
               if(out_size == out_max){
                  out_max *=2;
                  out_group=realloc(out_group,
                     out_max*sizeof(char*));
               }
               out_group[out_size++]=argv[optind++];
            }while(optind < argc && argv[optind][0]!='-');
            break;
         //'-g' option to specify which species to output an annotated genome
         //for, parsed in the same manner as -i and -o
         case 't':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "--tree-species parameter requires at least one argument\n");
               exit(1);
            }
            do{
               if(strcasestr(argv[optind],".maf")!=NULL) return;
               if(tree_species_size == tree_species_max){
                  tree_species_max *=2;
                  tree_species_names=realloc(tree_species_names,
                     tree_species_max*sizeof(char*));
               }
               tree_species_names[tree_species_size++]=argv[optind++];
            }while(optind < argc && argv[optind][0]!='-');
            break;
         //'-x' option specifies the in group conservation threshold
         case 'x':
            if(optarg==NULL || optarg[0]=='-'){
               fprintf(stderr, "--in-thresh parameter requires one argument\n");
               exit(1);
            }
            in_cons_thresh=atof(optarg);
            if(in_cons_thresh<=0 || in_cons_thresh >1){
               fprintf(stderr, "Invalid conservation threshold: %g\n",in_cons_thresh);
               exit(1);
            }
            break;
         //'-z' option specifies the out group conservation threshold
         case 'z':
            if(optarg==NULL || optarg[0]=='-'){
               fprintf(stderr, "--out-thresh parameter requires one argument\n");
               exit(1);
            }
            out_cons_thresh=atof(optarg);
            if(out_cons_thresh<=0 || out_cons_thresh >1){
               fprintf(stderr, "Invalid conservation threshold: %g\n",out_cons_thresh);
               exit(1);
            }
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
         //Default case to catch invalid options
         case '?':
	    if(optopt == NULL) fprintf(stderr,"Invalid long option: %s\n",argv[optind-1]);
	    else fprintf(stderr, "Invalid short option: %c\n", optopt);
            exit(1);
      }
   }
}

char **process_block_for_tree(sorted_alignment_block aln, char **tree_species, int num_species){
   //First initialize a ton of variables used throughout the loop, 'counts'
   //is an array of integers, which will store the count of each base seen 
   //in a given alignment column (6 options to allow for 'N' and '-')
   int counts[5] = {0};
   int itor = 0;
   int j = 0;
   int num_found = 0;
   int in_base = 0;
   int out_base = 0;
   int num_ones = 0;
   int found = 0;
   double in_score = 0.0;
   double out_score = 0.0;
   //'c' will hold the individual base being looked at in each column,
   //and cons_string will store the string of 0s, 1s, and 2s, as it is
   //built up one position at a time.
   char c;
   char *curr_species;
   int num_tree_species = 0;
   int num_seqs[num_species];
   for(;itor < num_species; ++itor) num_seqs[itor] = 0;
   for(itor = 0; itor < aln->in_size; ++itor){
      curr_species = (aln->in_sequences[itor])->species;
      for(;j < num_species; ++j)
         if(strcmp(curr_species,tree_species[j]) == 0){
            ++num_tree_species;  
            if(++num_seqs[j] == 2) return NULL;          
         }
   }
   for(itor = 0; itor < aln->out_size; ++itor){
      curr_species = (aln->out_sequences[itor])->species;
      for(j = 0;j < num_species; ++j)
         if(strcmp(curr_species,tree_species[j]) == 0){
            ++num_tree_species;
            if(++num_seqs[j] == 2) return NULL;
         }
   }
   if(num_tree_species < num_species - 2) return NULL;
   char **multiple_aln = malloc(sizeof(char *) * num_species);
   assert(multiple_aln != NULL);
   for(itor = 0; itor < num_species; ++itor){
      multiple_aln[itor] = calloc(aln->seq_length,sizeof(char));
      assert(multiple_aln[itor] != NULL); 
   }
   //Check conservation base by base, starting with in group species.
   for(unsigned int base = 0; base < aln->seq_length; ++base){
      //Need to reset the various variables used, particularly settings counts to 0
      itor=0;
      num_found=0;
      memset(counts,0,sizeof(counts));
      in_score=0.0;
      out_score=0.0;
      //Iterate over each in-group species in a given alignment column,
      //keeping counts of each base, note that N's aren't counted.
      for(;itor < aln->in_size;++itor){
         c=toupper(aln->in_sequences[itor]->sequence[base]);
         switch(c){
            case 'A':
                     ++counts[0];
                     break;
            case 'G':
                     ++counts[1];
                     break;
            case 'C':
                     ++counts[2];
                     break;
            case 'T':
                     ++counts[3];
                     break;
            case '-':
                     ++counts[4];
                     break;
            case 'N':
                     continue;
            default:
                     fprintf(stderr,"Nonstandard base encountered: %c\n",c);
              //return;
                     break;
         }
         ++num_found;
      }
//Get highest count found in this position, check if highest count over
//number of observed bases is below threshold, if so, continue, leaving
//the already written 0 untouched.
      if(num_found < 1){
          continue;
      }
      in_base=get_largest_index(counts,5);
      in_score=((double)counts[in_base])/num_found;
      if(in_score < in_cons_thresh){
         continue;
      }
//If in_score passes threshold, then check conservation in out group.
      itor = 0;
      num_found = 0;
      memset(counts,0,sizeof(counts));
      for(itor=0;itor < aln->out_size;++itor){
         c=toupper(aln->out_sequences[itor]->sequence[base]);
         switch(c){
            case 'A':
                     ++counts[0];
                     break;
            case 'G':
                     ++counts[1];
                     break;
            case 'C':
                     ++counts[2];
                     break;
            case 'T':
                     ++counts[3];
                     break;
            case '-':
                     ++counts[4];
                     break;
            case 'N':
                     continue;
            default:
              //fprintf(stderr,"Nonstandard base encountered: %c\n",c);
              //return;
              break;
         }
         ++num_found;
      }
      out_base=get_largest_index(counts,5);
      out_score=((double)counts[out_base])/num_found;
      if(num_found < 1 || in_base != out_base || out_score < out_cons_thresh){
         for(j=0;j < num_species ;++j){
            found = 0;
            curr_species = tree_species[j];
            for(int k = 0; k < aln->in_size; ++k)
               if(strcmp(curr_species, aln->in_sequences[k]->species) == 0){
                  multiple_aln[j][num_ones] = aln->in_sequences[k]->sequence[base];
                  found = 1;
                  break;
               }
            for(int k = 0; k < aln->out_size; ++k)
               if(strcmp(curr_species,aln->out_sequences[k]->species) == 0){
                  multiple_aln[j][num_ones] = aln->out_sequences[k]->sequence[base];
                  found = 1;
                  break;
               }
            if(!found) multiple_aln[j][num_ones] = 'N';
         }
         ++num_ones;
         continue;
      }
      //else cons_string[base]='2';
   }
   return multiple_aln;
}


void set_up(){
   in_group = malloc(sizeof(*in_group)*2);
   assert(in_group != NULL);
   in_size=0;
   in_max=2;
   out_group =malloc(sizeof(*out_group)*2);
   assert(out_group != NULL);
   out_size=0;
   out_max=2;
   in_cons_thresh=0.7;
   out_cons_thresh=0.7;
   tree_species_names = malloc(sizeof(char *)*2);
   tree_species_size=0;
   tree_species_max=2;
   num_columns = 0;

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
      fprintf(stderr, "Missing required MAF filename\n");
      exit(1);
   }
   char *filename = argv[optind];
   FILE *maf_file;
   if((maf_file= fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s",
         filename,strerror(errno));
      return 1;
   }
   int i = 0;
   for(; i < in_size; ++i) fprintf(stderr,"%s\n", in_group[i]);
   for(i = 0; i < out_size; ++i) fprintf(stderr,"%s\n", out_group[i]);
   for(i = 0; i < tree_species_size; ++i) fprintf(stderr,"%s\n", tree_species_names[i]);
   fprintf(stderr,"In Group Threshold: %g\n", in_cons_thresh);
   fprintf(stderr,"Out Group Threshold: %g\n", out_cons_thresh);
   fprintf(stderr,"Number of alignment columns: %d\n", num_columns);
   fprintf(stderr,"Filename: %s\n",filename);
   fprintf(stderr,"Num tree species: %d\n",tree_species_size);
   maf_linear_parser parser = get_linear_parser(maf_file,filename);
   char **multiple_alignment = malloc(sizeof(char *) * tree_species_size);
   assert(multiple_alignment != NULL);
   int curr_column = 0;
   int num_to_copy = 0;
   int aln_length = 0;
   for(i = 0; i < tree_species_size; ++i){
      multiple_alignment[i] = calloc(num_columns,sizeof(char));
      assert(multiple_alignment[i] != NULL);
   }
   while(1){
      sorted_alignment_block aln = get_sorted_alignment(parser,in_group
                  ,in_size,out_group,out_size);
      if(aln==NULL) break;
      char **aln_columns = process_block_for_tree(aln,tree_species_names,tree_species_size);
      if(aln_columns == NULL){
         free(aln_columns);
         free_sorted_alignment(aln);
         continue;
      }
      aln_length = strlen(aln_columns[0]);
      fprintf(stderr,"Aln length: %u\n",aln_length);
      if(aln_length == 0){ 
         free(aln_columns);
         free_sorted_alignment(aln);
         fprintf(stderr,"length zero\n");
         continue;
      }
      curr_column += aln_length;
      if(curr_column < num_columns) num_to_copy = aln->seq_length;
      else num_to_copy = (num_columns - curr_column);
      for(i = 0; i < tree_species_size; ++i){
         fprintf(stderr,"%u %s\n",strlen(aln_columns[i]),aln_columns[i]);
         strncat(multiple_alignment[i],aln_columns[i],num_to_copy);
         free(aln_columns[i]);
      }
      free(aln_columns);
      free_sorted_alignment(aln);
      fprintf(stderr,"One alignment block done");
      if(curr_column >= num_columns) break;
   }
   fprintf(stderr,"Num columns: %d\nCurr columns: %d\n",num_columns,curr_column);
   print_multiple_aln(multiple_alignment); 
   for(i = 0; i < tree_species_size; ++i)
      free(multiple_alignment[i]);
   free(multiple_alignment);
   free_linear_parser(parser);
   fclose(maf_file);
   clean_up();
   return 0;
}
