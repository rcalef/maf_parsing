#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <unistd.h>

#include "mafparser.h"

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
     clock_t start,end;
     double time_spent;
     start = clock();
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
        
     FILE *maf_file;
     if((maf_file= fopen(filename, "rb")) == NULL){
             fprintf(stderr, "Unable to open file: %s\nError: %s",
                            filename,strerror(errno));
             return 1;
     }
     maf_linear_parser parser = get_linear_parser(maf_file,filename);
//      while(1){
     alignment_block aln = linear_next_alignment_buffer(parser);
//           if(aln==NULL)break;
     print_alignment(aln);
     free_alignment_block(aln);
//        }
     aln = linear_next_alignment_buffer(parser);
     print_alignment(aln);
     free_alignment_block(aln);

     sorted_alignment_block saln=get_sorted_alignment(parser,
             in_group,in_size,out_group,out_size);
     print_sorted_alignment(saln);
     free_sorted_alignment(saln);
     free_linear_parser(parser);
     free(in_group);
     free(out_group);
     fclose(maf_file);
     end=clock();
     time_spent=(double)(end-start)/CLOCKS_PER_SEC;
     printf("Time spent: %g\n",time_spent);
     return 0;
}

