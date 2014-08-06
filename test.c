#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "mafparser.h"

int main(int argc, char **argv){
        clock_t start,end;
        double time_spent;
        start=clock();
        char *filename = argv[1];
        FILE *maf_file;
        if((maf_file= fopen(filename, "rb")) == NULL){
                fprintf(stderr, "Unable to open file: %s\nError: %s",
                                filename,strerror(errno));
                return 1;
        }
        maf_linear_parser parser = get_linear_parser(maf_file,filename);
        while(1){
           alignment_block aln = linear_next_alignment(parser);
           if(aln==NULL)break;
           print_alignment(aln);
           free_alignment_block(aln);
        }
        free_linear_parser(parser);
        fclose(maf_file);
        end=clock();
        time_spent=(double)(end-start)/CLOCKS_PER_SEC;
        printf("Time spent: %g\n",time_spent);
        return 0;
}

