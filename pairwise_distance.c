#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <time.h>


#include "mafparser.h"

typedef struct _pairwise_distance{
	int length;
	int num_idents;
	double percent;
}*dist;

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

dist get_pairwise_distance(seq seq1, seq seq2){
	unsigned int seq_length=strlen(seq1->sequence);
	if(seq_length != strlen(seq2->sequence)){
		fprintf(stderr,"Aligned sequences of different length\n"
			"Sequence 1: %s\nSequence 2: %s\n",
			seq1->sequence,seq2->sequence);
		return NULL;
	}
	dist distance=malloc(sizeof(*distance));
	distance->length = seq_length;
	distance->num_idents=0;
	for(unsigned int i=0;i<seq_length;++i)
		if(seq1->sequence[i]==seq2->sequence[i])
			++distance->num_idents;
	distance->percent = ((double)distance->num_idents)/distance->length;
	return distance;
}

void print_pairwise_distances(alignment_block aln){
	if(aln==NULL) return;
	dist distance;
	for(int i = 0; i < aln->size-1; ++i){
		for(int j=i+1;j<aln->size;++j){
			distance=get_pairwise_distance(aln->sequences[i],aln->sequences[j]);
			printf("Pairwise distance between %s and %s:\n"
				"%d/%d = %g\n",aln->sequences[i]->src,
				aln->sequences[j]->src,distance->num_idents,
				distance->length,distance->percent);
			free(distance);
		}
	}
}

int main(int argc, char **argv){
        clock_t start,end;
        double time_spent;
        start = clock();
        char *filename = argv[1];
        FILE *maf_file;
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
        printf("Time spent: %g\n",time_spent);
        return 0;
}

