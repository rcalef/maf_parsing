#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <search.h>

#include "mafparser.h"

typedef struct hsearch_data *hash;

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
   for(int i = 0; i < in_size; ++i) printf("%s\n", in_group[i]);
   for(int i = 0; i < out_size; ++i) printf("%s\n", out_group[i]);
   for(int i = 0; i < genomes_size; ++i) printf("%s\n", genome_names[i]);
   printf("Threshold: %g\n", cons_thresh);
   printf("Filename: %s\n",filename);
   genomes = calloc(1,sizeof(struct hsearch_data));
   int hc = hcreate_r(16,genomes);
   if(hc == 0){
      fprintf(stderr,"Failed to create hash table: %s\n", strerror(errno));
      exit(1);
   }ENTRY *ret_val;
   for(int i = 0; i < genomes_size; ++i){
//       ENTRY new_ent = malloc(sizeof(ENTRY));
//       new_ent.key=genome_names[i];
//       new_ent.data="This is the hashed data for " + genome_names[i]+"\n";
       char *data= malloc(sizeof(char) * 60);
       strcpy(data,"This is the hashed data for ");
       strcat(data,genome_names[i]);
       ENTRY new_ent = {genome_names[i],data};

       hc = hsearch_r(new_ent,ENTER,&ret_val,genomes);
       if(hc == 0){
          fprintf(stderr,"Failed to insert into hash table: %s\n", strerror(errno));
          exit(1);
        }
       printf("Entry inserted: %s\n", genome_names[i]);
   }
   for(int i = 0; i < genomes_size; ++i){
      ENTRY search={genome_names[i],NULL};
      hc = hsearch_r(search,FIND,&ret_val,genomes);
      if(hc == 0){
         fprintf(stderr,"Failed to read hash table: %s\n", strerror(errno));
         exit(1);
      }
      if(ret_val != NULL) printf("Entry found.\nKey: %s\nValue: %s\n",ret_val->key,ret_val->data);
   }
  
	   
/*        FILE *maf_file;
        if((maf_file= fopen(filename, "rb")) == NULL){
                fprintf(stderr, "Unable to open file: %s\nError: %s",
                                filename,strerror(errno));
                return 1;
        }
   int num_species= argc-2;
        char *species[num_species];
        for(int i =2; i < argc; ++i) species[i-2]=argv[i];
        maf_linear_parser parser = get_linear_parser(maf_file,filename);
        while(1){
           alignment_block aln = linear_next_alignment_buffer(parser);
           if(aln==NULL)break;
//           print_alignment(aln);
           species_filter(aln,species,num_species);
           print_alignment(aln);
      print_pairwise_distances(aln);
           free_alignment_block(aln);
        }
        free_linear_parser(parser);
        fclose(maf_file);
        end=clock();
        time_spent=(double)(end-start)/CLOCKS_PER_SEC;
        printf("Time spent: %g\n",time_spent);*/
        return 0;
}

