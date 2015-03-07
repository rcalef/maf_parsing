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

//The genome struct will store information about each genome that we want to
//output an annotated "genome" (fasta file of zeroes and ones) for. Each 
//struct stores a hash table of scaffolds, hashed by scaffold name.
//'num_scaffolds' and 'max_scaffolds' serve to keep track of the number of
//elements in the 'scaffolds' hash table, 'scaffold_names' serves to keep
//track of all scaffold names (and hence hash table keys), allowing iteration
//over the elements of the 'scaffolds' table. 'species' stores the name of
//the genome's species.
typedef struct _genome{
   int num_scaffolds;
   int max_scaffolds;
   hash scaffolds;
   char *species;
   char **scaffold_names;
}*genome;


//The scaffold struct simply stores the length of the sequence, and the
//sequence itself.
typedef struct _scaffold{
   unsigned int length;
   char *sequence;
}*scaffold;

typedef struct _thresh{
   double in_thresh;
   double out_thresh;
}*thresholds;


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
char **genome_names;
int genomes_size;
int genomes_max;
hash genomes;

//Define long options, note that options with 'no_argument'
//specified do require arguments, the 'no_argument' specification
//is necessary for reading in a variable size list of arguments.
static struct option long_options[]={
{"in-thresh",required_argument,0,'x'},
{"out-thresh",required_argument,0,'z'},
{"out-group",no_argument,0,'o'},
{"in-group",no_argument,0,'i'},
{"output-genomes",no_argument,0,'g'},
//Ignore the line below, necessary to define the long options
{0,0,0,0}
  };


//Convenience function to destroy a genome struct and free the memory used
void free_genome(genome gen){
   if(gen == NULL) return;
   
   free(gen->species);
   free(gen);
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
   free(in_group);
   free(out_group);
   ENTRY *gen_val;
   scaffold curr_scaf;
   genome curr_gen;
   for(int i = 0; i < genomes_size; ++i){
      gen_val=search_hash(genome_names[i],gen_val,genomes);
      curr_gen = gen_val->data;
      for(int j = 0; j < curr_gen->num_scaffolds; ++j){
         ENTRY *scaf_val;
         scaf_val = search_hash(curr_gen->scaffold_names[j],
             scaf_val,curr_gen->scaffolds);
         curr_scaf = scaf_val->data;
         free(curr_scaf->sequence);
         free(curr_scaf);
         free(scaf_val->key);
         free(curr_gen->scaffold_names[j]);
      }
      free(curr_gen->scaffold_names);
      hdestroy_r(curr_gen->scaffolds);
      free(curr_gen->scaffolds);
      free(curr_gen);
   }
   hdestroy_r(genomes);
   free(genomes);
   free(genome_names);
}

/*
species_filter is a function written to take in an alignment block, and filter
it out to only contain species of interest, not actually used in the body of
the program, but here in case it seems useful later
*/
void species_filter(alignment_block aln, char **species, 
             int num_species){
   if(aln==NULL) return;
   seq *new_sequences = malloc(sizeof(**new_sequences)*num_species);
   assert(new_sequences!=NULL);
   int size=0;
   char *species_name;
   char *temp;
   int cmp;
   for(int i=0; i < num_species;++i){
      for(int j = 0; j < aln->size;++j){
         temp=strdup(aln->sequences[j]->src);
         assert(temp!=NULL);
         species_name=strtok(temp,".");
         cmp=strcmp(species_name,species[i]);
         free(temp);
         if(cmp==0){
            new_sequences[size++]=copy_sequence(aln->sequences[j]);
            break;
         }
      }
   }
   for(int i =0; i < aln->size; ++i){
      free_sequence(aln->sequences[i]);
   }
   free(aln->sequences);
   aln->sequences=new_sequences;
   aln->size = size;
   aln->max = num_species;
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
   while((c=getopt_long(argc,argv,"x:z:iog",long_options,&option_index))!= -1){
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
         case 'g':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "--output-genomes parameter requires at least one argument\n");
               exit(1);
            }
            do{
               if(strcasestr(argv[optind],".maf")!=NULL) return;
               if(genomes_size == genomes_max){
                  genomes_max *=2;
                  genome_names=realloc(genome_names,
                     genomes_max*sizeof(char*));
               }
               genome_names[genomes_size++]=argv[optind++];
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
         //Default case to catch invalid options
         case '?':
	    if(optopt == NULL) fprintf(stderr,"Invalid long option: %s\n",argv[optind-1]);
	    else fprintf(stderr, "Invalid short option: %c\n", optopt);
            exit(1);
      }
   }
}

void check_and_initialize_scaffold(char *species, char *scaf_name, unsigned int scaf_size){
   if(!in_list(species,genome_names,
          genomes_size)) return;
   //If we are outputting this genome, then need to get the appropriate
   //genome struct.
   ENTRY *ret_val=search_hash(species,ret_val,genomes);
   genome curr_gen = ret_val->data;
   //Check if scaffold is in species genome struct already, if not, then
   //we need to add a new scaffold struct.
   ret_val=search_hash(scaf_name,ret_val,curr_gen->scaffolds);
   if(ret_val == NULL){
      //Must print warning if hash table begins to fill up, as the
      //C library hashing functions don't allow resizing the hash table.
      if(curr_gen->num_scaffolds >= curr_gen->max_scaffolds){
         fprintf(stderr, "WARNING: Scaffold hash table over half full"
                         " consider increasing max alignment hash size"
                         " to avoid decreased performance or crashes.\n"
                         "Species: %s\nCurrent size: %d\nMax size: %d\n"
                         ,curr_gen->species,curr_gen->num_scaffolds
                         ,curr_gen->max_scaffolds);
      }
      scaffold new_scaf= malloc(sizeof(*new_scaf));
      assert(new_scaf != NULL);
      new_scaf->length = scaf_size;
      new_scaf->sequence =  malloc(new_scaf->length*sizeof(char));
      assert(new_scaf->sequence != NULL);
      //When adding in a new scaffold, need to set the whole sequence to
      //0, corresponding to no information at unseen positions. Actually
      //needs to be the character 0, hence memset'ing to 48
      memset(new_scaf->sequence,48,new_scaf->length*sizeof(char));
      ENTRY search={strdup(scaf_name),new_scaf};
      assert(search.key != NULL);
      int hc=hsearch_r(search,ENTER,&ret_val,curr_gen->scaffolds);
      if(hc == 0){
         fprintf(stderr,"Error inserting into hash table: %s\n", strerror(errno));
         exit(1);
      }
      curr_gen->scaffold_names[curr_gen->num_scaffolds++]=
            strdup(scaf_name);
      assert(curr_gen->scaffold_names[curr_gen->num_scaffolds-1] != NULL);
   }
}


/*


*/
thresholds check_uninformative_block(sorted_alignment_block aln){
   unsigned int itor = 0;
   unsigned int num_in_species = 0;
   unsigned int num_out_species = 0;
   float in_species_proportion = 0.0;
   float out_species_proportion = 0.0;
   float new_in_thresh = 0.0;
   float new_out_thresh = 0.0;
   ENTRY *ret_val;
   hash in_group_species = calloc(1,sizeof(struct hsearch_data));
   assert(in_group_species != NULL);
   int hc = hcreate_r(4*in_size,in_group_species);
   if(hc == 0){
      fprintf(stderr,"Failed to create hash table: %s\n", strerror(errno));
      exit(1);
   }
   hash out_group_species = calloc(1,sizeof(struct hsearch_data));
   assert(out_group_species != NULL);
   hc = hcreate_r(4*out_size,out_group_species);
   if(hc == 0){
      fprintf(stderr,"Failed to create hash table: %s\n", strerror(errno));
      exit(1);
   }


   for(;itor < aln->in_size; ++itor){
      ret_val = search_hash(aln->in_sequences[itor]->species,
         ret_val,in_group_species);
      check_and_initialize_scaffold(aln->in_sequences[itor]->species, 
         aln->in_sequences[itor]->scaffold, aln->in_sequences[itor]->srcSize);
      if(ret_val != NULL) continue;
      ENTRY search = {aln->in_sequences[itor]->species,NULL};
      hc = hsearch_r(search,ENTER,&ret_val,in_group_species);
      if(hc == 0){
         fprintf(stderr,"Error inserting into hash table: %s\n", strerror(errno));
         exit(1);
      }
      ++num_in_species;
   }

   for(itor = 0;itor < aln->out_size; ++itor){
      ret_val = search_hash(aln->out_sequences[itor]->species,
         ret_val,out_group_species);
      if(ret_val != NULL) continue;
      ENTRY search = {aln->out_sequences[itor]->species,NULL};
      hc = hsearch_r(search,ENTER,&ret_val,out_group_species);
      if(hc == 0){
         fprintf(stderr,"Error inserting into hash table: %s\n", strerror(errno));
         exit(1);
      }
      ++num_out_species;
   }

   hdestroy_r(in_group_species);
   hdestroy_r(out_group_species);
   free(in_group_species);
   free(out_group_species);

   in_species_proportion = num_in_species/(float)in_size;
   printf("Num in species seen: %d\n In group species: %d\nIn group proportion: %g\n",
      num_in_species,in_size,in_species_proportion);
   if(in_species_proportion < in_cons_thresh) return NULL;
   new_in_thresh = (1 + in_cons_thresh) - in_species_proportion;


   out_species_proportion = num_out_species/(float)out_size;
   printf("Num out species seen: %d\n Out group species: %d\nOut group proportion: %g\n",
      num_out_species,out_size,out_species_proportion);
   //If out_species are 
   if(1 - out_species_proportion >= out_cons_thresh) new_out_thresh = -1;
   else if(out_species_proportion < out_cons_thresh) new_out_thresh = 2;
   new_out_thresh = (1 + out_cons_thresh) - out_species_proportion;


   thresholds new_thresh = malloc(sizeof(*new_thresh));
   assert(new_thresh != NULL);
   new_thresh->in_thresh = new_in_thresh;
   new_thresh->out_thresh = new_out_thresh;
   return new_thresh;
}

/*
process_block function is the meat of the program, this function actually
goes through each multiple alignment block, determining what string of
0s, 1s, and 2s, to output in the corresponding region of the genomes of 
interest.

Only input is a single 'sorted' alignment block, as defined in mafparser.h,
that is the alignment block object is pre-sorted into two sets of sequence
entries, one set for in-group species sequences, and one set for out-group
species sequences (sequences from species not in either group are thrown out).
*/
void process_block(sorted_alignment_block aln, thresholds thresh){
   //First initialize a ton of variables used throughout the loop, 'counts'
   //is an array of integers, which will store the count of each base seen 
   //in a given alignment column (6 options to allow for 'N' and '-')
   int counts[5] = {0};
   int itor = 0;
   int hc = 0;
   int num_found = 0;
   int offset = 0;
   int in_base = 0;
   int out_base = 0;
   double in_score = 0.0;
   double out_score = 0.0;
   ENTRY *ret_val;
   //'c' will hold the individual base being looked at in each column,
   //and cons_string will store the string of 0s, 1s, and 2s, as it is
   //built up one position at a time.
   char c;
   char cons_string[aln->seq_length];
   //thresh->in_thresh = 0.8;
   printf("Thresholds: %g %g\n",thresh->in_thresh,thresh->out_thresh);
   //Check conservation base by base, starting with in group species.
   for(unsigned int base = 0; base < aln->seq_length; ++base){
      //Need to reset the various variables used, particularly setting counts to 0
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
      //Check if we saw enough species to be convinced we have good information.
      //As of right now, cutoff implemented as 1, so right now serves to ignore
      //alignment blocks with no in-group species, which shouldn't happen
      //given that the alignment is archosaur referenced.

      //Get highest count found in this position, check if highest count over
      //number of observed bases is below threshold, if so, continue, leaving
      //the already written 0 untouched.
      if(num_found < 1){
          printf("Not enough found\n");
          cons_string[base]='0';
          continue;
      }
//********************************************************************
//Potential for rounding error here
//********************************************************************
      in_base=get_largest_index(counts,5);
      in_score=((float)counts[in_base])/num_found;
      printf("%g\n",in_score);
      if(in_score < thresh->in_thresh){
         printf("In score below thresh\n");
         cons_string[base]='0'; 
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
	      fprintf(stderr,"Nonstandard base encountered: %c\n",c);
	      //return;
	      break;
         }
         ++num_found;
      }

      //Like above, check if we have enough information on out-group species
      //if not, call it as in-group conserved (only get this far if the base
      //is in-group conserved).
      if(num_found < 1){
         cons_string[base]='1';
         continue;
      }
   
      //If the conserved base in the out-group is different from the conserved
      //base in the in-group, mark as in-group conserved. If bases are the same
      //then check if out-group conservation is above threshold, and mark base
      //as 1 if not, or 2 if so (conserved in both in-group and out-group)
      out_base=get_largest_index(counts,5);
      if(in_base != out_base) cons_string[base]='1';
      else{
         out_score=((float)counts[out_base])/num_found;
         if(out_score < thresh->out_thresh) cons_string[base]='1';
         else cons_string[base]='2';
      }
   }  

   //Now that we have the completed conservation string, we can add it
   //to the appropriate scaffold in the corresponding genome.
   //As of right now, have the unnecessary assumption that we're only 
   //outputting annotated genomes for in-group species, so we iterate
   //over the in-group species, continuing if the species is not in
   //the list of genomes to output.
   for(itor=0; itor < aln->in_size; ++itor){
      if(!in_list(aln->in_sequences[itor]->species,genome_names,
             genomes_size)) continue;
      //If we are outputting this genome, then need to get the appropriate
      //genome struct.
      ret_val=search_hash(aln->in_sequences[itor]->species,ret_val,genomes);
      genome curr_gen = ret_val->data;
      //Check if scaffold is in species genome struct already, if not, then
      //we need to add a new scaffold struct.
      ret_val=search_hash(aln->in_sequences[itor]->scaffold,ret_val,curr_gen->scaffolds);
/*
      if(ret_val == NULL){
         //Must print warning if hash table begins to fill up, as the
         //C library hashing functions don't allow resizing the hash table.
         if(curr_gen->num_scaffolds >= curr_gen->max_scaffolds){
            fprintf(stderr, "WARNING: Scaffold hash table over half full"
                            " consider increasing max alignment hash size"
                            " to avoid decreased performance or crashes.\n"
                            "Species: %s\nCurrent size: %d\nMax size: %d\n"
                            ,curr_gen->species,curr_gen->num_scaffolds
                            ,curr_gen->max_scaffolds);
         }
         scaffold new_scaf= malloc(sizeof(*new_scaf));
         assert(new_scaf != NULL);
         new_scaf->length = aln->in_sequences[itor]->srcSize;
         new_scaf->sequence =  malloc(new_scaf->length*sizeof(char));
         assert(new_scaf->sequence != NULL);
         //When adding in a new scaffold, need to set the whole sequence to
         //0, corresponding to no information at unseen positions. Actually
         //needs to be the character 0, hence memset'ing to 48
         memset(new_scaf->sequence,48,new_scaf->length*sizeof(char));
         ENTRY search={strdup(aln->in_sequences[itor]->scaffold),new_scaf};
	 assert(search.key != NULL);
         hc=hsearch_r(search,ENTER,&ret_val,curr_gen->scaffolds);
         if(hc == 0){
            fprintf(stderr,"Error inserting into hash table: %s\n", strerror(errno));
            exit(1);
         }
         curr_gen->scaffold_names[curr_gen->num_scaffolds++]=
               strdup(aln->in_sequences[itor]->scaffold);
         assert(curr_gen->scaffold_names[curr_gen->num_scaffolds-1] != NULL);
      }
*/
      //If scaffold entry already present, or after inserting new entry,
      //write to scaffold stream in appropriate position.
      
      //Only want to copy over the whole conservation string if the aligned
      //sequence for this species doesn't contain gaps, else need to only
      //copy over those numbers that correspond to existing bases.
      unsigned int insert_pos = aln->in_sequences[itor]->start;
      offset=0;
      //If the size of the sequence for this species is the same as the longest
      //sequence length, then just copy over the whole conservation string.
//*****************************************************************************
      //CHECK FOR CORRECTNESS 
//*****************************************************************************
      if(aln->in_sequences[itor]->size == aln->seq_length)
            memcpy(((scaffold)ret_val->data)->sequence+insert_pos,
                   cons_string,aln->seq_length*sizeof(char));
      //If sizes differ, go over the sequence for this species, skipping an
      //element of the conservation string whenever a '-' occurs in the sequence.
      else for(unsigned int i = 0; i < aln->seq_length; ++i){
	  if(aln->in_sequences[itor]->sequence[i] != '-'){
	     memcpy(((scaffold)ret_val->data)->sequence+insert_pos+offset,
                   cons_string+i,sizeof(char));
	     ++offset;
	  }
      }
   }
}

/*
write_genomes writes each genome struct to a separate fasta-like file,
line wrapping at 100 characters
*/
void write_genomes(){
   ENTRY *ret_val;
   char filename[64];
   FILE *outfile;
   for(int i = 0; i < genomes_size; ++i){
      strcpy(filename,genome_names[i]);
      strcat(filename,"_conservomatic.fasta");
      if((outfile= fopen(filename, "w")) == NULL){
         fprintf(stderr, "Unable to open file: %s\nError: %s",
            filename,strerror(errno));
         exit(1);
      }
      ret_val=search_hash(genome_names[i], ret_val,genomes);
      genome curr_gen = ret_val->data;
      for(int j = 0; j < curr_gen->num_scaffolds; ++j){
           ret_val=search_hash(curr_gen->scaffold_names[j],ret_val,curr_gen->scaffolds);
           fprintf(outfile,">%s.%s   ",genome_names[i],ret_val->key);
           char *sequence = ((scaffold)ret_val->data)->sequence;
           for(unsigned int k = 0; k < ((scaffold)ret_val->data)->length; ++k){
              if(k%100==0)fprintf(outfile,"\n");
              fprintf(outfile,"%c",sequence[k]);
           }
           fprintf(outfile,"\n");
      }
      fclose(outfile);
   }
}

void write_genomes_bed(){
   ENTRY *ret_val;
   char filename[64];
   FILE *outfile;
   int i = 0;
   int j = 0;
   unsigned int k = 0;
   int in_ones = 0;
   for(; i < genomes_size; ++i){
      strcpy(filename,genome_names[i]);
      strcat(filename,"_conservomatic.bed");
      if((outfile= fopen(filename, "w")) == NULL){
         fprintf(stderr, "Unable to open file: %s\nError: %s",
            filename,strerror(errno));
         exit(1);
      }
      ret_val=search_hash(genome_names[i], ret_val,genomes);
      genome curr_gen = ret_val->data;
      for(j = 0; j < curr_gen->num_scaffolds; ++j){
           ret_val=search_hash(curr_gen->scaffold_names[j],ret_val,curr_gen->scaffolds);
//           fprintf(outfile,">%s.%s   ",genome_names[i],ret_val->key);
           char *sequence = ((scaffold)ret_val->data)->sequence;
           for(k = 0; k < ((scaffold)ret_val->data)->length; ++k){
              if(sequence[k] == '1'){
                 if(in_ones) continue;
                 else{
                    fprintf(outfile, "%s\t%u\t", ret_val->key,k);
                    in_ones = 1;
                    continue;
                 }
              }
              else if(in_ones){
                    fprintf(outfile,"%u\n",k);
                    in_ones = 0;
                    continue;
              }
           }
           if(in_ones){
              fprintf(outfile,"%u\n",k);
              in_ones = 0;
           }
      }
      fclose(outfile);
   }

}


/*
print_genomes is just like write_genomes, execpt printing to stdout
*/
void print_genomes(){
   ENTRY *ret_val;
   for(int i = 0; i < genomes_size; ++i){
      printf("For species %s:\n",genome_names[i]);
      ret_val=search_hash(genome_names[i], ret_val,genomes);
      genome curr_gen = ret_val->data;
      for(int j = 0; j < curr_gen->num_scaffolds; ++j){
           ret_val=search_hash(curr_gen->scaffold_names[j],ret_val,curr_gen->scaffolds);
           printf(">%s.%s   ",genome_names[i],ret_val->key);
           char *sequence = ((scaffold)ret_val->data)->sequence;
           for(unsigned int k = 0; k < ((scaffold)ret_val->data)->length; ++k){
              if(k%100==0)printf("\n");
              printf("%c",sequence[k]);
           }
           printf("\n");
      }
   }
}


//Set up is a simple convenience function to initialize the global
//variables to default values.
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
   genome_names = malloc(sizeof(char *)*2);
   genomes_size=0;
   genomes_max=2;
   genomes = calloc(1,sizeof(struct hsearch_data));
   assert(genomes != NULL);
   int hc = hcreate_r(16,genomes);
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
   if((maf_file= fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s",
         filename,strerror(errno));
      return 1;
   }

   //Print out command line arguments as a sanity check
   for(int i = 0; i < in_size; ++i) printf("%s\n", in_group[i]);
   for(int i = 0; i < out_size; ++i) printf("%s\n", out_group[i]);
   for(int i = 0; i < genomes_size; ++i) printf("%s\n", genome_names[i]);
   printf("In Group Threshold: %g\n", in_cons_thresh);
   printf("Out Group Threshold: %g\n", out_cons_thresh);
   printf("Filename: %s\n",filename);
   printf("In group size: %d\n",in_size);
   printf("Out group size: %d\n",out_size); 


   //Next we need to initialize the genome structs for each genome we
   //want to output. Assume a large number of scaffolds as we can't
   //resize the hash tables after creation.
   ENTRY *ret_val;
   int hc = 0;
   for(int i = 0; i < genomes_size; ++i){
       genome new_gen = malloc(sizeof(*new_gen));
       assert(new_gen != NULL);
       new_gen->num_scaffolds=0;
       new_gen->max_scaffolds=600000;
       //Set max to less than half the actual max, as we don't want the
       //hash table more than half full.
       new_gen->scaffold_names = malloc(sizeof(char*) * 2000000);
       assert(new_gen->scaffold_names != NULL);
       new_gen->species = genome_names[i];
       new_gen->scaffolds = calloc(1,sizeof(struct hsearch_data));
       assert(new_gen->scaffolds != NULL);
       hc = hcreate_r(2000000,new_gen->scaffolds);
       if(hc == 0){
          fprintf(stderr,"Failed to create hash table: %s\n", strerror(errno));
          exit(1);
       }
    
       ENTRY new_ent = {genome_names[i],new_gen};
       hc = hsearch_r(new_ent,ENTER,&ret_val,genomes);
       if(hc == 0){
          fprintf(stderr,"Failed to insert into hash table: %s\n", strerror(errno));
          exit(1);
        }
       printf("Entry inserted: %s\n", genome_names[i]);
   }

   //After genome structs are all set up, just need to get the MAF parser object
   //and then process each alignmnent block, one at a time.
   thresholds thresh = NULL;
   maf_linear_parser parser = get_linear_parser(maf_file,filename);
   while(1){
      sorted_alignment_block aln = get_sorted_alignment(parser,in_group
                  ,in_size,out_group,out_size);
      //Parser returns NULL on EOF
      if(aln==NULL)break;
      thresh = check_uninformative_block(aln);
      if(thresh == NULL){
         free_sorted_alignment(aln);
         continue;
      }
      process_block(aln,thresh);
      free_sorted_alignment(aln);
   }
   //Finally just write genomes to files and clean up memory usage.
   write_genomes();
   write_genomes_bed();
   free_linear_parser(parser);
   fclose(maf_file);
   clean_up();
   return 0;
}

