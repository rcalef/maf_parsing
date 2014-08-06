#ifndef __MAFPARSER_H
#define __MAFPARSER_H

#define BUFSIZE 50000
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
}*seq;

typedef struct _alignment_block{
	double score;
	int pass;
	char *data;
	seq *sequences;
	int size;
	int max;
        int curr_seq;
}*alignment_block;

int get_next_offset(maf_array_parser parser);
seq get_sequence(char *data);
seq copy_sequence(seq sequence);

alignment_block array_next_alignment(maf_array_parser parser);
alignment_block linear_next_alignment(maf_linear_parser parser);
alignment_block linear_next_alignment_buffer(maf_linear_parser parser);
void array_double(maf_array_parser parser);

maf_array_parser get_array_parser(FILE *maf_file,char *filename);
maf_linear_parser get_linear_parser(FILE *maf_file, char *filename);
void free_array_parser(maf_array_parser parser);
void free_linear_parser(maf_linear_parser parser);
void free_sequence(seq sequence);
void free_alignment_block(alignment_block aln);

void print_sequence(seq sequence);
void print_alignment(alignment_block aln);
#endif
