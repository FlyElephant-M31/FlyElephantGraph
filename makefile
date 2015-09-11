# makefile for programs and generator

# compiler flags:
#
# CHECK_SOLUTION     check if at the end the sink is unreachable 
#                    from the source
# CUT_ONLY           define to compute min-cut but not max-flow
# PRINT_STAT         define to print detailed operation statistics
# PRINT_FLOW         define to print out flow value
# WAVE_INIT          wave-style initialization (initial global update,
#                    then a wave of flow pushing). Better on some problems
# OLD_INIT           define to use old-style initialization (large excess
#                    at the source, global update)
# INIT_UPDATE        perform initial global update
# EXCESS_TYPE_LONG   set excessType to long, otherwise it is long long int
#                    if "long long int" not supported, change to "double"
#                    (in types.h)
#
# VERBOSE            print info on gus: s, t, thread numbers
# DONTPERMUTE        do not randomly permute the nodes     
# DEGREE_HEURISTIC   order nodes by degree
# DECREASING	     order nodes by non-increasing degree
# DEGREE_HEURISTIC   check if d(s) = flow. If so, do not change the tree
# NOTRIVIALCUT       dont make the trivial cut test

CCOMP = gcc
GUSFLAGS = -O3 -Wall -DSHORTFLOW -DDONTPERMUTE -DDEGREE -DDEGREE_HEURISTIC -DDECREASING 

all: gus_hipr gusp

gus_hipr: gus_hipr.c parser_undirected.c types.h timer.c
	$(CCOMP) $(GUSFLAGS) -o gus_hipr gus_hipr.c

# gusp.c is gusp_opt3.c
gusp: gusp.c parser_undirected.c types.h timer.c
	$(CCOMP) $(GUSFLAGS) -fopenmp -o gusp gusp.c	

clean: 
	rm -f gus_hipr gusp *.o *~
