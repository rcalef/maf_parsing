#define _GNU_SOURCE
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <ctype.h>

int main(int argc, char **argv){
   char *in_filename = argv[1];
   FILE *bed_file;
   if((bed_file= fopen(in_filename, "r")) == NULL){
      fprintf(stderr, "Unable to open file: %s\nError: %s",
         in_filename,strerror(errno));
      return 1;
   }
   char buffer[1024];
   char geneName[256];
   char *nl;
   unsigned int start=0;
   unsigned int end =0;
   unsigned long total=0;
   char buf[1024];
   while(fscanf(bed_file,"%s\t%u\t%u\t%s",buffer,&start,&end,buf) ==4){
//Remove newline from buffer
      if((nl=strchr(buffer,'\n'))!=NULL)
         *nl = '\0';
      total += (end - start);
   }
   printf("Total number of bases covered: %lu\n",total);
   fclose(bed_file);
   return 0;
}

