#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <string.h>
#include <assert.h>
#include <time.h>

#include "mafparser.h"

int main(int argc, char **argv){
	int *nums = calloc(10,sizeof(int));
	for(int i = 0; i < 10; ++i) printf("%d",nums[i]);
	printf("\n");
        return 0;
}

