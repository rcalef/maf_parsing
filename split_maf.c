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

//First define structs that will be used to represent various data structures


//Next define global variables used to keep track of parameters passed in to
char *split_species;
char **scaffold_names;
int names_size;
int names_max;
int scaffolds_size;
int scaffolds_max;
hash scaffold_outputs;

//Define long options, note that options with 'no_argument'
//specified do require arguments, the 'no_argument' specification
//is necessary for reading in a variable size list of arguments.
static struct option long_options[]={
{"split-species",required_argument,0,'s'},
//Ignore the line below, necessary to define the long options
{0,0,0,0}
  };




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
   ENTRY *gen_val;
   for(int i = 0; i < scaffolds_size; ++i){
      gen_val=search_hash(scaffold_names[i],gen_val,scaffold_outputs);
      free(gen_val->data);
   }
   for(int i = 0; i < names_size; ++i) free(scaffold_names[i]);
   hdestroy_r(scaffold_outputs);
   free(scaffold_names);
   free(scaffold_outputs);
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
   while((c=getopt_long(argc,argv,"s",long_options,&option_index))!= -1){
      switch(c){
         //'-i' option to specify in group species, takes in an arbitrary sized
         //list of species
         case 's':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "--species parameter requires exactly one argument\n");
               exit(1);
            }
            if(strcasestr(argv[optind],".maf")!=NULL){
               fprintf(stderr,"--species paramater requires exactly one argument\n");
               exit(1);
            }
               
            split_species=argv[optind++];
            break;
         //Default case to catch invalid options
         case '?':
	    if(optopt == NULL) fprintf(stderr,"Invalid long option: %s\n",argv[optind-1]);
	    else fprintf(stderr, "Invalid short option: %c\n", optopt);
            exit(1);
      }
   }
}


void process_block(alignment_block aln){
   ENTRY *ret_val;
   seq split_species_seq = NULL;
   char filename[128];
   char *curr_scaffold = NULL;
   char *scaf_filename = NULL;
   int found = 0;
   int hc = 0;
   for(int i = 0; i < aln->size; ++i){
      if(strcmp(aln->sequences[i]->species,split_species) == 0){
         split_species_seq = aln->sequences[i];
         found = 1;
         break;
      }
   }
   if(!found) curr_scaffold = "not_present";
   else curr_scaffold = split_species_seq->scaffold;
   ret_val = search_hash(curr_scaffold,ret_val,scaffold_outputs);
   if(ret_val == NULL){
      if(scaffolds_size >= scaffolds_max){
         fprintf(stderr, "WARNING: Scaffold hash table over half full"
                         " consider increasing max alignment hash size"
                         " to avoid decreased performance or crashes.\n"
                         "Current size: %d\nMax size: %d\n"
                         ,scaffolds_size,scaffolds_max);
       }
       strcpy(filename,split_species);
       strcat(filename,".");
       strcat(filename,curr_scaffold);
       strcat(filename,"-blocks_only.maf");
       curr_scaffold = strdup(curr_scaffold);
       assert(curr_scaffold != NULL);
       ENTRY search = {curr_scaffold,strdup(filename)};
       assert(search.key != NULL);
       hc=hsearch_r(search,ENTER,&ret_val,scaffold_outputs);
       if(hc == 0){
          fprintf(stderr,"Error inserting into hash table: %s\n", strerror(errno));
          exit(1);
       }
       if(names_size >= names_max){
          names_max *= 2;
          scaffold_names = realloc(scaffold_names,names_max * sizeof(char *));
       }
       ++scaffolds_size; 
       scaffold_names[names_size++] = curr_scaffold;
   }
   scaf_filename = ret_val->data;
   FILE *curr_scaf_output = fopen(scaf_filename,"a");
   if(curr_scaf_output == NULL){
      fprintf(stderr,"Unable to open file: %s\nError: %s\n",
         filename,strerror(errno));
      exit(1);
   }
   print_alignment_file(curr_scaf_output,aln);
   fclose(curr_scaf_output);
}

//Set up is a simple convenience function to initialize the global
//variables to default values.
void set_up(){
   split_species = NULL;
   names_size = 0;
   names_max = 2048;
   scaffold_names = malloc(sizeof(char *) * names_max);
   assert(scaffold_names != NULL);
   scaffold_outputs = calloc(1,sizeof(struct hsearch_data));
   assert(scaffold_outputs != NULL);
   scaffolds_size = 0;
   scaffolds_max = 8192;
   int hc = hcreate_r(4*scaffolds_max,scaffold_outputs);
   if(hc == 0){
      fprintf(stderr,"Failed to create hash table: %s\n", strerror(errno));
      exit(1);
   }
}

int main(int argc, char **argv){
   //First we need to set up the global alignments, and parse the command
   //line arguments, opening the appropriate file for reading.
   set_up();
   parse_args(argc,argv);
   if(optind >= argc){
      fprintf(stderr, "Missing required MAF filename\n");
      exit(1);
   }
   char *filename = argv[optind];
   FILE *maf_file;
   if((maf_file = fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s",
         filename,strerror(errno));
      return 1;
   }
   maf_linear_parser parser = get_linear_parser(maf_file,filename);
   while(1){
      alignment_block aln = linear_next_alignment_buffer(parser);
      if(aln == NULL) break;
      //Parser returns NULL on EOF
      process_block(aln);
      free_alignment_block(aln);
   }
   //Finally just write genomes to files and clean up memory usage.
   free_linear_parser(parser);
   fclose(maf_file);
   clean_up();
   return 0;
}

