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

typedef struct _genome{
   int num_scaffolds;
   int max_scaffolds;
   hash scaffolds;
   char *species;
   char **scaffold_names;
}*genome;

typedef struct _scaffold{
   unsigned int length;
   char *sequence;
}*scaffold;


char **in_group;
int in_size;
int in_max;
char **out_group;
int out_size;
int out_max;
double in_cons_thresh;
double out_cons_thresh;
char **genome_names;
int genomes_size;
int genomes_max;
hash genomes;

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
 



void free_genome(genome gen){
   if(gen == NULL) return;
   
   free(gen->species);
   free(gen);
}

int get_largest(int *nums, int size){
   int max = 0;
   for(int i =0; i < size; ++i) 
      if(nums[i] > max) max = nums[i];
   return max;
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

void clean_up(){
   free(in_group);
   free(out_group);
   ENTRY *gen_val;
//   ENTRY *scaf_val;
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

void parse_args(int argc, char **argv){
   if (argc < 4){
// print_usage();
      fprintf(stderr,"Too few arguments\n");
      exit(1);
   }
   char c;
   int option_index=0;
   while((c=getopt_long(argc,argv,"x:z:iog",long_options,&option_index))!= -1){
      switch(c){
         case 'i':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "--in-group parameter requires at least one argument\n");
               exit(1);
            }
            do{
               if(strcasestr(argv[optind],".maf")!=NULL) return;
               if(in_size == in_max){
                  in_max*=2;
                  in_group=realloc(in_group,
                     in_max*sizeof(char*));
               }
               in_group[in_size++]=argv[optind++];
            }while(optind < argc && argv[optind][0]!='-');
            break;
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
         case '?':
	   //	   if(optopt == NULL) fprintf(stderr,"Invalid long option: %s\n",argv[optind-1]);
	   //           else fprintf(stderr, "Invalid short option: %s\n", optopt);
            exit(1);
      }
   }
}


void process_block(sorted_alignment_block aln){
   int counts[5] = {0};
   int itor;
   int hc = 0;
   int num_found;
   int offset;
   double in_score;
   double out_score;
   ENTRY *ret_val;
   char c;
   char cons_string[aln->seq_length];
//   memset(cons_string,0,aln->seq_length*sizeof(int));
//Check conservation base by base, starting with in group species.
   for(unsigned int base = 0; base < aln->seq_length; ++base){
      itor=0;
      num_found=0;
      memset(counts,0,sizeof(counts));
      in_score=0.0;
      out_score=0.0;
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
	      //fprintf(stderr,"Nonstandard base encountered: %c\n",c);
  	      //return;
                     break;        
         }                   
         ++num_found;
      }  
//Get highest count found in this position, check if highest count over
//number of observed bases is below threshold, if so, continue, leaving
//the already written 0 untouched.
      if(num_found < 1){
          cons_string[base]='0';
          continue;
      }
      in_score=((double)get_largest(counts,5))/num_found;
      if(in_score < in_cons_thresh){
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
	      //fprintf(stderr,"Nonstandard base encountered: %c\n",c);
	      //return;
	      break;
         }
         ++num_found;
      }
      if(num_found < 1){
         cons_string[base]='1';
         continue;
      }
      out_score=((double)get_largest(counts,5))/num_found;
      if(out_score < out_cons_thresh) cons_string[base]='1';
      else cons_string[base]='2';
   }  
//Now that we have the completed conservation string, we can add it
//to the appropriate scaffold in the corresponding genome.
   for(itor=0; itor < aln->in_size; ++itor){
//First check if species genome is being outputted.
//If not, continue.
      if(!in_list(aln->in_sequences[itor]->species,genome_names,
             genomes_size)) continue;
//If so, get scaffold name and genome struct.
      ret_val=search_hash(aln->in_sequences[itor]->species,ret_val,genomes);
      genome curr_gen = ret_val->data;
//Check if scaffold is in species genome struct already.
      ret_val=search_hash(aln->in_sequences[itor]->scaffold,ret_val,curr_gen->scaffolds);
//If ret_val is NULL, need to add entry for this scaffold
      if(ret_val == NULL){
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
//If scaffold entry already present, or after inserting new entry,
//write to scaffold stream in appropriate position.
//      print_sequence(aln->in_sequences[itor]);
      unsigned int insert_pos = aln->in_sequences[itor]->start;
//Only want to copy over the whole conservation string if the aligned
//sequence for this species doesn't contain gaps, else need to only
//copy over those numbers that correspond to existing bases.
      offset=0;
      if(aln->in_sequences[itor]->size == aln->seq_length)
            memcpy(((scaffold)ret_val->data)->sequence+insert_pos,
                   cons_string,aln->seq_length*sizeof(char));
      else for(int i = 0; i < aln->seq_length; ++i){
	  if(aln->in_sequences[itor]->sequence[i] != '-'){
	     memcpy(((scaffold)ret_val->data)->sequence+insert_pos+offset,
                   cons_string+i,sizeof(char));
	     ++offset;
	  }
      }
   }
}

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

int main(int argc, char **argv){
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
   for(int i = 0; i < in_size; ++i) printf("%s\n", in_group[i]);
   for(int i = 0; i < out_size; ++i) printf("%s\n", out_group[i]);
   for(int i = 0; i < genomes_size; ++i) printf("%s\n", genome_names[i]);
   printf("In Group Threshold: %g\n", in_cons_thresh);
   printf("Out Group Threshold: %g\n", out_cons_thresh);
   printf("Filename: %s\n",filename);
   genomes = calloc(1,sizeof(struct hsearch_data));
   assert(genomes != NULL);
   int hc = hcreate_r(16,genomes);
   if(hc == 0){
      fprintf(stderr,"Failed to create hash table: %s\n", strerror(errno));
      exit(1);
   }
   ENTRY *ret_val;
   for(int i = 0; i < genomes_size; ++i){
       genome new_gen = malloc(sizeof(*new_gen));
       assert(new_gen != NULL);
       new_gen->num_scaffolds=0;
       new_gen->max_scaffolds=600000;
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
   maf_linear_parser parser = get_linear_parser(maf_file,filename);
   while(1){
      sorted_alignment_block aln = get_sorted_alignment(parser,in_group
                  ,in_size,out_group,out_size);
      if(aln==NULL)break;
      process_block(aln);
      free_sorted_alignment(aln);
   }
  // print_genomes();
   write_genomes();
/*   for(int i = 0; i < genomes_size; ++i){
      printf("For species %s:\n",genome_names[i]);
      ret_val=search_hash(genome_names[i], ret_val,genomes);
      genome curr_gen = ret_val->data;
      for(int j = 0; j < curr_gen->num_scaffolds; ++j){
           ret_val=search_hash(curr_gen->scaffold_names[j],ret_val,curr_gen->scaffolds);
           printf("\t%s   %s\n",ret_val->key,(char *)ret_val->data);
      }
   }*/
   free_linear_parser(parser);
   fclose(maf_file);
   clean_up();
   return 0;
}

