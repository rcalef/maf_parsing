/*
 * maf_parser.c
 *
 *  Created on: Aug 2, 2014
 *      Author: calef_000
 */

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>

struct parser{
	int alignment_blocks[10];
	int curr_block;
        int size;
};

typedef struct parser *maf_parse;

int get_next_offset(maf_parse parser){
   if(parser->curr_block >= parser->size) return -1;
   return parser->alignment_blocks[parser->curr_block++];
}

int main(int argc, char **argv){
	char *filename = argv[1];
	FILE *maf_file;
	if((maf_file= fopen(filename, "rb")) == NULL){
		fprintf(stderr, "Unable to open file: %s\nError: %s",
				filename,strerror(errno));
		return 1;
	}
	char buffer[3000];
	assert(buffer != NULL);
	maf_parse parser = (maf_parse *)malloc(sizeof(maf_parse));
	assert(parser != NULL);
	int pos;
	int i =-1;
	while(!feof(maf_file)){
		pos = ftell(maf_file);
		char *check = fgets(buffer,3000,maf_file);
		if(ferror(maf_file) != 0){
			fprintf(stderr, "File stream error: %s\nError: %s",
					filename,strerror(errno));
			return 1;
		}
		if(buffer[0]=='a'){
			parser->alignment_blocks[++i]=pos;
		}
	}
        parser->size=i+1;
        parser->curr_block=0;
	int j=0;
	for(; j <= i; ++j) printf("%d\n",parser->alignment_blocks[j]);
        j=0;
        while(1){
           int offset = get_next_offset(parser);
           if(offset==-1)break;
           int check =fseek(maf_file, offset,SEEK_SET);
           if(check!=0)fprintf(stderr, "File seek error: %s\nError: %s",
                         filename,strerror(errno));
           char *fc = fgets(buffer,3000,maf_file);
           printf("Alignment block %d\n%s\n",++j,buffer);
        }
	return 0;
}
