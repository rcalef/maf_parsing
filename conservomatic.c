#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <ctype.h>

#include "mafparser.h"

typedef struct _genome{
   int num_scaffolds;
   int max_scaffolds;
   hash scaffolds;
   char *species;
}*genome;


char **in_group;
int in_size;
int in_max;
char **out_group;
int out_size;
int out_max;
double cons_thresh;
char **genome_names;
int genomes_size;
int genomes_max;
hash genomes;

int get_largest(int *nums, int size){
   int max = 0;
   for(int i =0; i < size; ++i) 
      if(nums[i] > max) max = nums[i];
   return max;
}

void search_hash(char *key, ENTRY *ret_val, hash table){
         ENTRY search={key,NULL};
         int hc = hsearch_r(search,FIND,&ret_val,table);
         if(hc == 0){
            fprintf(stderr,"Error searching hash table: %s\n", strerror(errno));
            exit(1);
         }
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
   while((c=getopt(argc,argv,"ciog"))!= -1){
      switch(c){
         case 'i':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "-i parameter requires at least one argument\n");
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
               fprintf(stderr, "-o parameter requires at least one argument\n");
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
               fprintf(stderr, "-g parameter requires at least one argument\n");
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
         case 'c':
            if(argv[optind][0]=='-'){
               fprintf(stderr, "-c parameter requires one argument\n");
               exit(1);
            }
            cons_thresh=atof(argv[optind++]);
            if(cons_thresh<=0 || cons_thresh >1){
               fprintf(stderr, "Invalid conservation threshold: %g\n",cons_thresh);
               exit(1);
            }
            break;
         case '?':
            fprintf(stderr, "Invalid option: %s\n", argv[optind-1]);
            exit(1);
      }
   }
}


void process_block(sorted_alignment_block aln){
   int counts[5] = {0};
   int itor =0;
   int hc = 0;
   int num_found = 0;
   double in_score = 0.0;
   double out_score = 0.0;
   ENTRY *ret_val;
   char c;
   char cons_string[aln->seq_length];
   memset(cons_string,0,aln->seq_length*sizeof(char));
//Check conservation base by base, starting with in group species.
   for(unsigned int base = 0; base < aln->seq_length; ++base){
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
            default:
                    fprintf(stderr,"Nonstandard base encountered: %c\n",c);
                    return;
         }                   
         ++num_found;
      }  
//Get highest count found in this position, check if highest count over
//number of observed bases is below threshold, if so, continue, leaving
//the already written 0 untouched.
      in_score=((double)get_largest(counts,5))/num_found;
      if(in_score < cons_thresh) continue; 
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
            default:
                    fprintf(stderr,"Nonstandard base encountered: %c\n",c);
                    return;
         }
         ++num_found;
      }
      out_score=((double)get_largest(counts,5))/num_found;
      if(out_score < cons_thresh) cons_string[base]=1;
      else cons_string[base]=2;
   }  
//Now that we have the completed conservation string, we can add it
//to the appropriate scaffold in the corresponding genome.
   char *scaffold_name;
   char *temp;
   for(itor=0; itor < genomes_size; ++itor){
//First check if species was in alignment block.
//      ENTRY search={genome_names[itor],NULL};
//      hc=hsearch_r(search,FIND,&ret_val,aln->sequences);
//      if(hc == 0){
//         fprintf(stderr,"Error searching hash table: %s\n", strerror(errno));
//         exit(1);
//      }
      search_hash(genome_names[itor],ret_val,aln->sequences);
      if(ret_val == NULL) continue;
      seq curr_seq = ret_val->data;
//If species was present, get scaffold name and genome struct.
//      hc=hsearch_r(search,FIND,&ret_val,genomes);
//      if(hc == 0){
//         fprintf(stderr,"Error searching hash table: %s\n", strerror(errno));
//         exit(1);
//      }
      search_hash(genome_names[itor],ret_val,genomes);
      genome curr_gen = ret_val->data;
      temp = strdup(((seq)ret_val->data)->src);
      scaffold_name=strtok(temp,".");
      scaffold_name=strtok(NULL,".");
      printf("Scaffold name: %s\n",scaffold_name);
//Check if scaffold is in species genome struct already.
//      search.key=scaffold_name;
//      hc=hsearch_r(search,FIND,&ret_val,curr_gen->scaffolds);
//      if(hc == 0){
//         fprintf(stderr,"Error searching hash table: %s\n", strerror(errno));
//         exit(1);
//      }
      search_hash(scaffold_name,ret_val,curr_gen->scaffolds);
//If ret_val is NULL, need to add entry for this scaffold
      if(ret_val == NULL){
         char *scaffold_stream =  calloc(1,sizeof(char)*curr_seq->srcSize);
         assert(scaffold_stream != NULL);
         ENTRY search={scaffold_name,scaffold_stream};
         hc=hsearch_r(search,ENTER,&ret_val,curr_gen->scaffolds);
         if(hc == 0){
            fprintf(stderr,"Error searching hash table: %s\n", strerror(errno));
            exit(1);
         }
      }
//If scaffold entry already present, or after inserting new entry,
//write to scaffold stream in appropriate position.
      unsigned int insert_pos = curr_seq->start;
      strcpy(((char*)ret_val->data)[insert_pos],cons_string);
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
   cons_thresh=0;
        parse_args(argc,argv);
   if(optind >= argc){
      fprintf(stderr, "Missing required MAF filename\n");
      exit(1);
   }
   char *filename = argv[optind];
//   for(int i = 0; i < in_size; ++i) printf("%s\n", in_group[i]);
//  for(int i = 0; i < out_size; ++i) printf("%s\n", out_group[i]);
//  for(int i = 0; i < genomes_size; ++i) printf("%s\n", genome_names[i]);
//   printf("Threshold: %g\n", cons_thresh);
//   printf("Filename: %s\n",filename);
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
       new_gen->max_scaffolds=128;
       new_gen->species = genome_names[i];
       new_gen->scaffolds = calloc(1,sizeof(struct hsearch_data));
       assert(new_gen->scaffolds != NULL);
       hc = hcreate_r(256,new_gen->scaffolds);
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
	   
   FILE *maf_file;
   if((maf_file= fopen(filename, "rb")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s",
         filename,strerror(errno));
      return 1;
   }
   maf_linear_parser parser = get_linear_parser(maf_file,filename);
   while(1){
      sorted_alignment_block aln = get_sorted_alignment(parser,in_group
                  ,in_size,out_group,out_size);
      if(aln==NULL)break;
           print_sorted_alignment(aln);
           process_block(aln);
           free_sorted_alignment(aln);
        }
        free_linear_parser(parser);
        fclose(maf_file);
        return 0;
}

