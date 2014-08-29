# $Id: Makefile,v 1.8 2014-05-14 18:07:06-07 - - $

MKFILE    = Makefile
NOINCLUDE = ci clean spotless
NEEDINCL  = ${filter ${NOINCLUDE}, ${MAKECMDGOALS}}

GCC       = gcc -g -O0 -Wall -Wextra -std=gnu99
MKDEPS    = gcc -MM

STATSSOURCE = maf_stats.c mafparser.c
STATSOBJECTS = ${STATSSOURCE: .c=.o}
CONSSOURCE   = conservomatic.c mafparser.c
CONSOBJECTS   = ${CSOURCE:.c=.o}
EXECBIN   = conservomatic maf_stats
SOURCES   = ${CHEADER} ${CSOURCE} ${MKFILE}
TESTCMD   = ./conservomatic --in-group amaVit1 croPor2 Anc05 Anc14 Anc21 --out-group Anc10 Anc09 Anc07 Anc18 --in-thresh=0.8 --out-thresh=0.7 --output-genomes amaVit1 croPor2 Anc05 larger_artificial.maf


all : ${EXECBIN}

conservomatic: ${CONSOBJECTS}
	${GCC} -o $@ ${CONSOBJECTS}

maf_stats : ${STATSOBJECTS}
	${GCC} -o $@ ${STATSOBJECTS}

%.o : %.c
	${GCC} -c $<


clean :
	- rm ${OBJECTS} ${DEPSFILE} core ${EXECBIN}.errs

test: ${EXECBIN}
	${TESTCMD} > test.out
	diff Anc05_conservomatic_testcheck.fasta Anc05_conservomatic.fasta > test.check
	diff amaVit1_conservomatic_testcheck.fasta amaVit1_conservomatic.fasta >> test.check
	diff croPor2_conservomatic_testcheck.fasta croPor2_conservomatic.fasta >> test.check
again :
	${GMAKE} spotless deps ci all lis

ifeq "${NEEDINCL}" ""
include ${DEPSFILE}
endif
