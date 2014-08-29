#ifndef __MAFPARSER_H
#define __MAFPARSER_H

#include <search.h>

#define BUFSIZE 50000

typedef struct hsearch_data *hash;

typedef struct array_parser{
        FILE *maf_file;
        char *filename;
        int *alignment_blocks;
        int curr_block;
        int size;
        int max;
}*maf_array_parser;

typedef struct linear_parser{
	FILE *maf_file;
	char *filename;
	char buf[BUFSIZE];
        unsigned long curr_pos;
        int fill_buf;
        char *pos;
}*maf_linear_parser;

typedef struct _aligned_sequence{
	char *src;
	unsigned long start;
	unsigned int size;
	char strand;
	unsigned long srcSize;
	char *sequence;
        char *species;
	char *scaffold;
}*seq;

typedef struct _alignment_block{
	double score;
	int pass;
	char *data;
	seq *sequences;
	int size;
	int max;
        int curr_seq;
	unsigned int seq_length;
}*alignment_block;

typedef struct _sorted_alignment_block{
	double score;
	int pass;
	char *data;
	unsigned int seq_length;
	int in_size;
	int in_max;
	seq *in_sequences;
	int out_size;
	int out_max;
	seq *out_sequences;
}*sorted_alignment_block;

typedef struct _hash_alignment_block{
        double score;
        int pass;
        char *data;
	int size;
        int max;
        unsigned int seq_length;
	char **species;
	hash sequences;
}*hash_alignment_block;


int in_list(char *needle, char **haystack, int size);
void array_double(maf_array_parser parser);

int get_next_offset(maf_array_parser parser);
seq get_sequence(char *data);
seq copy_sequence(seq sequence);

alignment_block array_next_alignment(maf_array_parser parser);
alignment_block linear_next_alignment(maf_linear_parser parser);
alignment_block linear_next_alignment_buffer(maf_linear_parser parser);
hash_alignment_block get_next_alignment_hash(maf_linear_parser parser);
sorted_alignment_block get_sorted_alignment(maf_linear_parser parser,
              char **in_group, int in_size, char **out_group, int out_size);

maf_array_parser get_array_parser(FILE *maf_file,char *filename);
maf_linear_parser get_linear_parser(FILE *maf_file, char *filename);
void free_array_parser(maf_array_parser parser);
void free_linear_parser(maf_linear_parser parser);
void free_sequence(seq sequence);
void free_alignment_block(alignment_block aln);
void free_sorted_alignment(sorted_alignment_block aln);
void free_hash_alignment(hash_alignment_block aln);

void print_sequence(seq sequence);
void print_alignment(alignment_block aln);
void print_sorted_alignment(sorted_alignment_block aln);
void print_hash_alignment(hash_alignment_block aln);
#endif
