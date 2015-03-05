#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <ctype.h>

#define THRESH 500
#define MIN_LENGTH 100

int main(int argc, char **argv){
   char *in_filename = argv[1];
   FILE *cons_file;
   if((cons_file= fopen(in_filename, "r")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s",
         in_filename,strerror(errno));
      return 1;
   }
   char *out_filename = argv[2];
   FILE *bed_file;
   if((bed_file= fopen(out_filename, "w")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s",
         out_filename,strerror(errno));
      return 1;
   }
   char buffer[1024];
   char *curr_chr=0;
   char *nl;
   unsigned int pos=0;
   unsigned int line_length=0;
   unsigned int start=0;
   unsigned int end=0;
   unsigned int break_count=0;
   unsigned int block_length=0;
   short in_block=0;
   while(fgets(buffer,1024,cons_file) != NULL){
      if(buffer[0]=='>'){
         if(in_block){
            end=pos-1;
            block_length -= break_count;
            if(block_length >= MIN_LENGTH)
               fprintf(bed_file,"%s,%u,%u\n",curr_chr,start,end);
            block_length=0;
            in_block=0;
            break_count=0;
         }
         free(curr_chr);
         curr_chr=strtok(buffer,".\n");
         curr_chr=strdup(strtok(NULL," .\n"));
         pos=0;
         continue;
      }
//Remove newline from buffer
      if((nl=strchr(buffer,'\n'))!=NULL)
         *nl = '\0';
      line_length=strlen(buffer);
      for(unsigned int i = 0; i < line_length; ++i){
         if(!in_block){
            if(buffer[i]=='1'){
               start=pos+i;
               in_block=1;
               continue;
            }
            else continue;
         }
//Only want to break if we see THRESH number of 0s or 2s consecutively
//so set break_count to 0 if we see a 1
         if(buffer[i] != '1') ++break_count;
         else break_count=0;
         ++block_length;
         if(break_count == THRESH){
            end = pos + i - break_count;
            block_length -= break_count;
            if(block_length >= MIN_LENGTH)
               fprintf(bed_file,"%s,%u,%u\n",curr_chr,start,end);
            block_length=0;
            in_block=0;
            break_count=0;
         }
     }
     pos += line_length;
   }
   free(curr_chr);
   fclose(cons_file);
   fclose(bed_file);
   return 0;
}

