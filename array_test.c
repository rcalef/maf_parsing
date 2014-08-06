#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

#include "mafparser.h"

int main(int argc, char **argv){
        char *filename = argv[1];
        FILE *maf_file;
        if((maf_file= fopen(filename, "rb")) == NULL){
                fprintf(stderr, "Unable to open file: %s\nError: %s",
                                filename,strerror(errno));
                return 1;
        }
        maf_array_parser parser = get_array_parser(maf_file,filename);
        int j=0;
        for(; j < parser->size; ++j) printf("%d\n",parser->alignment_blocks[j]);
        j=0;
        char buffer[3000];
        while(1){
           int offset = get_next_offset(parser);
           if(offset==-1)break;
           int check =fseek(maf_file, offset,SEEK_SET);
           if(check!=0)fprintf(stderr, "File seek error: %s\nError: %s",
                         filename,strerror(errno));
           char *fc = fgets(buffer,3000,maf_file);
           printf("Alignment block %d\n%s\n",++j,buffer);
        }
        parser->curr_block=1;
        int offset=get_next_offset(parser);
           int check =fseek(maf_file, offset,SEEK_SET);
           if(check!=0)fprintf(stderr, "File seek error: %s\nError: %s",
                         filename,strerror(errno));
          char *fc = fgets(buffer,3000,maf_file);
          fc = fgets(buffer,3000,maf_file);
        printf("%s\n",buffer);
        alignment_block aln = get_next_alignment(parser);
//        seq test=aln->sequences[0];
//        printf("Alignment 1 sequence 1: \n%s  %lu  %u  %c  %lu  %s\n"
//            ,test->src,test->start,test->size,test->strand,
//                test->srcSize,test->sequence);
        free_alignment_block(aln);
        free_parser(parser);
        fclose(maf_file);
        return 0;
}

