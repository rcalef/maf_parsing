#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <ctype.h>
#include <getopt.h>

#include "mafparser.h"

typedef struct _block{
   unsigned int num_blocks;
   unsigned int *sequence_counts;
   unsigned int num_counts;
   unsigned int max_counts;
   unsigned int *species_counts;
   unsigned int num_species;
   unsigned int max_species;
}*block_stats;

typedef struct _species{
   char *species;
   unsigned int *seqs_per_block;
   unsigned int num_seqs;
   unsigned int max_seqs;
   unsigned int *length_per_block;
   unsigned int num_lengths;
   unsigned int max_lengths;
}*species_stats;

hash total_species_stats;
char **species_in_stats;
int num_spec;
block_stats block;
hash temp_counts;
char **species_seen;
unsigned int num_species_seen;

//Define long options, note that options with 'no_argument'
//specified do require arguments, no_argument specification
//necessary for reading in variable size list of arguments.
static struct option long_options[]={
{"in-thresh",required_argument,0,'x'},
{"out-thresh",required_argument,0,'z'},
{"out-group",no_argument,0,'o'},
{"in-group",no_argument,0,'i'},
{"output-genomes",no_argument,0,'g'},
{0,0,0,0}
  };
 
block_stats new_block_stats(){
   block_stats block = malloc(sizeof(struct _block));
   assert(block != NULL);
   block->num_blocks = 0;
   block->sequence_counts = malloc(sizeof(unsigned int) * 256);
   assert(block->sequence_counts != NULL);
   block->num_counts=0;
   block->max_counts=256;
   block->species_counts = malloc(sizeof(unsigned int) * 256);
   assert(block->species_counts != NULL);
   block->num_species=0;
   block->max_species=256;
   return block;
}

species_stats new_species_stats(char *species_name){
   species_stats species = malloc(sizeof(struct _species));
   assert(species != NULL);
   species->species = strdup(species_name);
   assert(species->species != NULL);
   species->seqs_per_block = malloc(sizeof(unsigned int) * 256);
   assert(species->seqs_per_block != NULL);
   species->num_seqs=0;
   species->max_seqs=256;
   species->length_per_block = malloc(sizeof(unsigned int) * 256);
   assert(species->length_per_block != NULL);
   species->num_lengths=0;
   species->max_lengths=256;
   return species;
}

ENTRY *search_hash(char *key, ENTRY *ret_val, hash table){
         ENTRY search={key,NULL};
         int hc = hsearch_r(search,FIND,&ret_val,table);
         if(hc == 0){
            if(errno == 3) return NULL;
            fprintf(stderr,"Error searching hash table: %s\n", strerror(errno));
            exit(1);
         }
         return ret_val;
}



void process_block(alignment_block aln){
   ENTRY *ret_val;
   species_stats curr_stats;
   seq curr_seq;
   unsigned int curr_count;
   num_species_seen=0;
   int hc = hcreate_r(256,temp_counts);
   if(hc == 0){
      fprintf(stderr,"Failed to create hash table: %s\n", strerror(errno));
      exit(1);
   }
   printf("%d\n",aln->size);
   for(int i = 0; i < aln->size; ++i){
      curr_seq = aln->sequences[i];
//Check if species already has an entry in counts table, if not, add one
      ret_val = search_hash(curr_seq->species,ret_val,temp_counts);
      if(ret_val == NULL){
          unsigned int *count;
          *count=0;
          ENTRY insert={curr_seq->species,count};
          hc = hsearch_r(insert,ENTER,&ret_val,temp_counts);
          if(hc == 0){
             fprintf(stderr,"Error inserting into hash table: %s\n", strerror(errno));
             exit(1);
          }
          species_seen[num_species_seen++]=curr_seq->species;
        //  printf("Entry inserted for species: %s\n", curr_seq->species);    
      }
      ++(*((unsigned int *)ret_val->data));
    //  printf("%u\n",(*((unsigned int *)ret_val->data)));
   }
//Next we need to add the temp counts to our overall counts.
   for(unsigned int i = 0; i < num_species_seen; ++i){
//First get the count.
      printf("%s\n",species_seen[i]);
      ret_val = search_hash(species_seen[i],ret_val,temp_counts);
      curr_count = (*(unsigned int *)ret_val->data);
//Check if species already has an entry in stats hash table, if not, add an entry
      ret_val = search_hash(species_seen[i],ret_val,total_species_stats);
      if(ret_val == NULL){
          species_stats stats = new_species_stats(curr_seq->species);
          ENTRY insert={strdup(species_seen[i]),stats};
          hc = hsearch_r(insert,ENTER,&ret_val,total_species_stats);
          if(hc == 0){
             fprintf(stderr,"Error inserting into hash table: %s\n", strerror(errno));
             exit(1);
          }
          species_in_stats[num_spec++]=strdup(species_seen[i]);
          printf("New species found: %s\n", species_in_stats[num_spec-1]);
      }
      else printf("Already in there: %s\n", species_seen[i]);
      curr_stats = ((species_stats)ret_val->data);
//Insert new count of sequences in block to end of array,
//doubling array if necessary
      if(curr_stats->num_seqs == curr_stats->max_seqs){
         curr_stats->max_seqs *= 2;
         curr_stats->seqs_per_block=realloc(curr_stats->seqs_per_block,
            curr_stats->max_seqs*sizeof(unsigned int));
     }
     curr_stats->seqs_per_block[curr_stats->num_seqs++] = curr_count;
//Similarly for length of sequences in the block
     if(curr_stats->num_lengths == curr_stats->max_lengths){
        curr_stats->max_lengths *= 2;
        curr_stats->length_per_block=realloc(curr_stats->length_per_block, 
           curr_stats->max_lengths*sizeof(unsigned int));
     }
     curr_stats->length_per_block[curr_stats->num_lengths++]=aln->seq_length;
   }
//Adjust block stats, and destroy current temp count table.
   ++block->num_blocks;
   if(block->num_counts == block->max_counts){
     block->max_counts *= 2;
     block->sequence_counts=realloc(block->sequence_counts,
        block->max_counts*sizeof(unsigned int));
   }
   block->sequence_counts[block->num_counts++] = aln->size;
   if(block->num_species == block->max_species){
      block->max_species *= 2;
      block->species_counts=realloc(block->species_counts,
         block->max_species*sizeof(unsigned int));
   }
   block->species_counts[block->num_species++] = num_species_seen;
   hdestroy_r(temp_counts);
}

double get_variance(unsigned int *values, 
       unsigned int num_values, double mean){
   double variance = 0;
   for(unsigned int i = 0 ; i < num_values; ++i)
      variance += ((values[i] - mean) * (values[i]-mean));
   return variance/num_values;
}

void print_block_stats(block_stats stats){
   unsigned int total_seqs=0;
   unsigned int total_species=0;
   printf("Number of blocks: %u\nNumber of sequences per block:\n", stats->num_blocks);
   for(unsigned int i = 0; i < stats->num_counts; ++i){
      printf("%u\n",stats->sequence_counts[i]);
      total_seqs += stats->sequence_counts[i];
   }
   printf("Number of species per block:\n");
   for(unsigned int i =0; i < stats->num_species; ++i){
      printf("%u\n",stats->species_counts[i]);
      total_species += stats->species_counts[i];
   }
   double seq_average= ((double)total_seqs)/stats->num_blocks;
   double species_average = ((double)total_species)/stats->num_blocks;
   printf("Average number of sequences per block: %g\n",seq_average);
   printf("   Variance: %g\n", get_variance(stats->sequence_counts,
      stats->num_counts,seq_average));
   printf("Average number of species per block: %g\n",species_average);
   printf("   Variance: %g\n", get_variance(stats->species_counts,
      stats->num_species,species_average));
}

void print_species_stats(species_stats stats){
   unsigned int total_seqs=0;
   unsigned int total_lengths=0;
   printf("For species %s\n   Number of sequences per block:\n",
      stats->species);
   for(unsigned int i = 0; i < stats->num_seqs; ++i){
      printf("   %u\n",stats->seqs_per_block[i]);
      total_seqs += stats->seqs_per_block[i];
   }
   printf("   Length per block:\n");
   for(unsigned int i =0; i < stats->num_lengths; ++i){
      printf("   %u\n",stats->length_per_block[i]);
      total_lengths += stats->length_per_block[i];
   }
   double seq_average=((double)total_seqs)/block->num_blocks;
   printf("   Average number of sequences per block: %g\n",seq_average);
   printf("   Variance: %g\n", 
      get_variance(stats->seqs_per_block,stats->num_seqs,seq_average));
   double length_average=((double)total_lengths)/stats->num_seqs;
   printf("   Average length of sequences: %g\n",length_average);
   printf("   Variance: %g\n",
      get_variance(stats->length_per_block,stats->num_lengths,length_average));
}

int main(int argc, char **argv){
   char *filename = argv[1];
   FILE *maf_file;
   if((maf_file= fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s",
         filename,strerror(errno));
      return 1;
   }
   printf("Filename: %s\n",filename);
   block = new_block_stats();
   total_species_stats = calloc(1, sizeof(struct hsearch_data));
   assert(total_species_stats != NULL);
   int hc = hcreate_r(256,total_species_stats);
   if(hc == 0){
      fprintf(stderr,"Failed to create hash table: %s\n", strerror(errno));
      exit(1);
   }
   temp_counts = calloc(1,sizeof(struct hsearch_data));
   assert(temp_counts != NULL);
   species_seen = calloc(100,sizeof(char *));
   num_species_seen = 0;
   species_in_stats=calloc(100,sizeof(char *));
   num_spec=0;
   maf_linear_parser parser = get_linear_parser(maf_file,filename);
   while(1){
      alignment_block aln = linear_next_alignment_buffer(parser);
      if(aln==NULL)break;
      process_block(aln);
      free_alignment_block(aln);
   }
   print_block_stats(block);
   ENTRY *ret_val;
   for(int i = 0; i < num_spec; ++i){
      ret_val=search_hash(species_in_stats[i],ret_val,total_species_stats);
      print_species_stats((species_stats)ret_val->data);
  //    printf("%s\n",species_in_stats[i]);
   }
   free_linear_parser(parser);
   fclose(maf_file);
   //clean_up();
   return 0;
}

