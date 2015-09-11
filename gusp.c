/* 
   This code implements a parallel version of Gusfield's flow equivalent tree 
   algorithm using OpenMP.

   If results are reported, references should include:

   J. Cohen, L. A. Rodrigues, F. Silva, R. Carmo, A. Guedes, E. P. Duarte Jr., 
   "Parallel Implementations of Gusfield's Cut Tree Algorithm," 
   11th International Conference Algorithms and Architectures for Parallel Processing (ICA3PP), 
   pp. 258-269, Lecture Notes in Computer Science (LNCS) 7016, ISSN 0302-9743, Melbourne,
   Australia, 2011. 

   This code is derived from HIPR by IG Systems, Inc.
   HIPR implements a maximum flow - highest level push-relabel algorithm 
   http://www.igsystems.com/hipr/
   Commercial use requires a license.
*/

/* gusp_opt3.c */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <time.h>
#include <omp.h> //(lar)

#include "values.h"
#include "types.h"                     /* type definitions */
#include "parser_undirected.c"         /* parser */
#include "timer.c"                     /* timing routine */

/*
#define GLOB_UPDT_FREQ 0.5
*/
#define GLOB_UPDT_FREQ 0.5
#define ALPHA 6
#define BETA 12

#define WHITE 0
#define GREY 1
#define BLACK 2

/* global variables */

long   n;                    /* number of nodes */
long   m;                    /* number of arcs */
long   nm;                   /* n + ALPHA * m */
long   nMin;                 /* smallest node id */
cType  *cap;                 /* array of capacities */

float globUpdtFreq;          /* global update frequency */
int seed;

// typedef struct graph *Graph;
typedef struct graph {
node   *nodes;               /* array of nodes */
arc    *arcs;                /* array of arcs */
bucket *buckets;             /* array of buckets */

node   *source;              /* source node pointer */
node   *sink;                /* sink node pointer */
//node   **queue;              /* queue for BFS */
//node   **qHead, **qTail, **qLast;     /* queue pointers */
long   dMax;                 /* maximum label */
long   aMax;                 /* maximum active node label */
long   aMin;                 /* minimum active node label */
double flow;                 /* flow value */
long pushCnt;           /* number of pushes */
long relabelCnt;       /* number of relabels */
long updateCnt;       /* number of updates */
long gapCnt;           /* number of gaps */
long gNodeCnt;           /* number of nodes after gap */  
float t, t2;                 /* for saving times */
node   *sentinelNode;        /* end of the node list marker */
arc *stopA;                  /* used in forAllArcs */
long workSinceUpdate;      /* the number of arc scans since last update */
  /* more variables that can't be shared by threads */
long i_dist;
node *i_next, *i_prev;

} *Graph;//end graph

// begin: variables declarations (jc)
long *tree; 
double *tree_weights;
long step; 
// end: variables (jc) 

/* macros */

#define forAllNodes(i, G) for ( i = G->nodes; i != G->sentinelNode; i++ )
#define forAllArcs(i,a, G) for (a = i->first, G->stopA = (i+1)->first; a != G->stopA; a++)

#define nNode( i, G ) ( (i) - G->nodes + nMin )
// Index(): from pointer to index (jc)
#define Index( i, G ) ( (i) - G->nodes )
#define nArc( a, G )  ( ( a == NULL )? -1 : (a) - G->arcs )

#define min( a, b ) ( ( (a) < (b) ) ? a : b )

/* FIFO queue for BFS macros */
/*
#define qInit() \
{\
  qHead = qTail = queue;\
}

#define qEmpty ( qHead == qTail )

#define qEnqueue(i) \
{\
  *qTail = i;\
  if ( qTail == qLast ) qTail = queue;\
  else qTail++;\
}

#define qDequeue(i) \
{\
  i = *qHead;\
  if ( qHead == qLast ) qHead = queue;\
  else qHead++;\
}
*/

/* 
   bucket macros:
   bucket's active node list is singly-linked
     operations aAdd, aRemove (from the front)
   bucket's inactive list is doubly-linked
     operations iAdd, iDelete (from arbitrary position)
*/


#define aAdd(l,i, G)\
{\
  i->bNext = l->firstActive;\
  l->firstActive = i;\
  G->i_dist = i->d;\
  if (G->i_dist < G->aMin)\
    G->aMin = G->i_dist;\
  if (G->i_dist > G->aMax)\
    G->aMax = G->i_dist;\
  if (G->dMax < G->aMax)\
    G->dMax = G->aMax;\
}

/* i must be the first element */
#define aRemove(l,i)				\
{\
  l->firstActive = i->bNext;\
}

#define iAdd(l,i, G)\
{\
  G->i_next = l->firstInactive;\
  i->bNext = G->i_next;\
  i->bPrev = G->sentinelNode;\
  G->i_next->bPrev = i;\
  l->firstInactive = i;\
}

#define iDelete(l,i, G)				\
{\
  G->i_next = i->bNext;\
  if (l->firstInactive == i) {\
    l->firstInactive = G->i_next;\
    G->i_next->bPrev = G->sentinelNode;\
  }\
  else {\
    G->i_prev = i->bPrev;\
    G->i_prev->bNext = G->i_next;\
    G->i_next->bPrev = G->i_prev;\
  }\
}

/* allocate datastructures, initialize related variables */

int allocDS( Graph G )
{

  nm = ALPHA * n + m;
  /*
  queue = (node**) calloc ( n, sizeof (node*) );
  if ( queue == NULL ) return ( 1 );
  qLast = queue + n - 1;
  qInit();
  */
  G->buckets = (bucket*) calloc ( n+2, sizeof (bucket) );
  if ( G->buckets == NULL ) return ( 1 );

  G->sentinelNode = G->nodes + n;
  G->sentinelNode->first = G->arcs + 2*m;

  return ( 0 );

} /* end of allocate */


void init( Graph G )
{
  node  *i;        /* current node */
  int overflowDetected;
  bucket *l;
  arc *a;
#ifdef EXCESS_TYPE_LONG
  double testExcess;
#endif
#ifndef OLD_INIT
  unsigned long delta;
#endif

   G->pushCnt  = 0;           /* number of pushes */
   G->relabelCnt   = 0;       /* number of relabels */  
   G->updateCnt    = 0;       /* number of updates */
   G->gapCnt   = 0;           /* number of gaps */
   G->gNodeCnt = 0;           /* number of nodes after gap */  
   G->workSinceUpdate=0;  

  // initialize excesses

  forAllNodes(i, G) {
    i->excess = 0;
    i->current = i->first;
    forAllArcs(i, a, G)
      {
	a->resCap = cap[a-G->arcs];
	// a->rev->resCap = cap[a-arcs]; // for undirected graphs, did not work. jc
      }
  }

  for (l = G->buckets; l <= G->buckets + n-1; l++) {
    l -> firstActive   = G->sentinelNode;
    l -> firstInactive  = G->sentinelNode;
  }
    
  overflowDetected = 0;
#ifdef EXCESS_TYPE_LONG
  testExcess = 0;
  forAllArcs(source,a, G) {
    if (a->head != G->source) {
      testExcess += a->resCap;
    }
  }
  if (testExcess > MAXLONG) {
    printf("c WARNING: excess overflow. See README for details.\nc\n");
    overflowDetected = 1;
  }
#endif
#ifdef OLD_INIT
  G->source -> excess = MAXLONG;
#else
  if (overflowDetected) {
    G->source -> excess = MAXLONG;
  }
  else {
    G->source->excess = 0;
    forAllArcs(G->source,a, G) {
      if (a->head != G->source) {
	G->pushCnt ++;
	delta = a -> resCap;
	a -> resCap -= delta;
	(a -> rev) -> resCap += delta;
	a->head->excess += delta;
      }
    }
  }

  /*  setup labels and buckets */
  l = G->buckets + 1;
    
  G->aMax = 0;
  G->aMin = n;
    
  forAllNodes(i, G) {
    if (i == G->sink) {
      i->d = 0;
      iAdd(G->buckets,i, G);
      continue;
    }
    if ((i == G->source) && (!overflowDetected)) {
      i->d = n;
    }
    else
      i->d = 1;
    if (i->excess > 0) {
      /* put into active list */
      aAdd(l,i, G);
    }
    else { /* i -> excess == 0 */
      /* put into inactive list */
      if (i->d < n)
	iAdd(l,i, G);
    }
  }
  G->dMax = 1;
#endif

  //  dMax = n-1;
  G->flow = 0.0;

} /* end of init */

void checkMax(Graph G)
{
  bucket *l;

  for (l = G->buckets + G->dMax + 1; l < G->buckets + n; l++) {
    assert(l->firstActive == G->sentinelNode);
    assert(l->firstInactive == G->sentinelNode);
  }
}

/* global update via backward breadth first search from the sink */

void globalUpdate (Graph G)
{

  node  *i, *j;       /* node pointers */
  arc   *a;           /* current arc pointers  */
  bucket *l, *jL;          /* bucket */
  long curDist, jD;
  long state;


  G->updateCnt ++;

  /* initialization */

  forAllNodes(i, G)
    i -> d = n;
  G->sink -> d = 0;

  // Era "l <= buckets + dmax".
  // Essa modificacao evita o segmentation fault na
  // instancia dblcyc.1024 e produz a resposta correta.
  // As linha abaixo apenas inicializam os buckets.
  for (l = G->buckets; l <= G->buckets + n; l++) {
    l -> firstActive   = G->sentinelNode;
    l -> firstInactive  = G->sentinelNode;
  }

  G->dMax = G->aMax = 0;
  G->aMin = n;

  /* breadth first search */

  // add sink to bucket zero

  iAdd(G->buckets, G->sink, G);
  for (curDist = 0; 1; curDist++) {

    state = 0;
    l = G->buckets + curDist;
    jD = curDist + 1;
    jL = l + 1;
    /*
    jL -> firstActive   = sentinelNode;
    jL -> firstInactive  = sentinelNode;
    */

    if ((l->firstActive == G->sentinelNode) && 
	(l->firstInactive == G->sentinelNode))
      break;

    while (1) {

      switch (state) {
      case 0: 
	i = l->firstInactive;
	state = 1;
	break;
      case 1:
	i = i->bNext;
	break;
      case 2:
	i = l->firstActive;
	state = 3;
	break;
      case 3:
	i = i->bNext;
	break;
      default: 
	assert(0);
	break;
      }
      
      if (i == G->sentinelNode) {
	if (state == 1) {
	  state = 2;
	  continue;
	}
	else {
	  assert(state == 3);
	  break;
	}
      }

      /* scanning arcs incident to node i */
      forAllArcs(i,a, G) {
	if (a->rev->resCap > 0 ) {
	  j = a->head;
	  if (j->d == n) {
	    j->d = jD;
	    j->current = j->first;
	    if (jD > G->dMax) G->dMax = jD;
	    
	    if (j->excess > 0) {
	      /* put into active list */
	      aAdd(jL,j, G);
	    }
	    else {
	      /* put into inactive list */
	      iAdd(jL,j, G);
	    }
	  }
	}
      } /* node i is scanned */ 
    }
  }

} /* end of global update */

/* second stage -- preflow to flow */
void stageTwo ( Graph G )
/*
   do dsf in the reverse flow graph from nodes with excess
   cancel cycles if found
   return excess flow in topological order
*/

/*
   i->d is used for dfs labels 
   i->bNext is used for topological order list
   buckets[i-nodes]->firstActive is used for DSF tree
*/

{
  node *i, *j, *tos, *bos, *restart, *r;
  arc *a;
  cType delta;

  /* deal with self-loops */
  forAllNodes(i, G) {
    forAllArcs(i,a, G)
      if ( a -> head == i ) {
	a -> resCap = cap[a - G->arcs];
      }
  }

  /* initialize */
  tos = bos = NULL;
  forAllNodes(i, G) {
    i -> d = WHITE;
    //    buckets[i-nodes].firstActive = NULL;
    G->buckets[i-G->nodes].firstActive = G->sentinelNode;
    i -> current = i -> first;
  }

  /* eliminate flow cycles, topologicaly order vertices */
  forAllNodes(i, G)
    if (( i -> d == WHITE ) && ( i -> excess > 0 ) &&
	( i != G->source ) && ( i != G->sink )) {
      r = i;
      r -> d = GREY;
      do {
	for ( ; i->current != (i+1)->first; i->current++) {
	  a = i -> current;
	  if (( cap[a - G->arcs] == 0 ) && ( a -> resCap > 0 )) { 
	    j = a -> head;
	    if ( j -> d == WHITE ) {
	      /* start scanning j */
	      j -> d = GREY;
	      G->buckets[j-G->nodes].firstActive = i;
	      i = j;
	      break;
	    }
	    else
	      if ( j -> d == GREY ) {
		/* find minimum flow on the cycle */
		delta = a -> resCap;
		while ( 1 ) {
		  delta = min ( delta, j -> current -> resCap );
		  if ( j == i )
		    break;
		  else
		    j = j -> current -> head;
		}

		/* remove delta flow units */
		j = i;
		while ( 1 ) {
		  a = j -> current;
		  a -> resCap -= delta;
		  a -> rev -> resCap += delta;
		  j = a -> head;
		  if ( j == i )
		    break;
		}
	  
		/* backup DFS to the first saturated arc */
		restart = i;
		for ( j = i -> current -> head; j != i; j = a -> head ) {
		  a = j -> current;
		  if (( j -> d == WHITE ) || ( a -> resCap == 0 )) {
		    j -> current -> head -> d = WHITE;
		    if ( j -> d != WHITE )
		      restart = j;
		  }
		}
	  
		if ( restart != i ) {
		  i = restart;
		  i->current++;
		  break;
		}
	      }
	  }
	}

	if (i->current == (i+1)->first) {
	  /* scan of i complete */
	  i -> d = BLACK;
	  if ( i != G->source ) {
	    if ( bos == NULL ) {
	      bos = i;
	      tos = i;
	    }
	    else {
	      i -> bNext = tos;
	      tos = i;
	    }
	  }

	  if ( i != r ) {
	    i = G->buckets[i-G->nodes].firstActive;
	    i->current++;
	  }
	  else
	    break;
	}
      } while ( 1 );
    }


  /* return excesses */
  /* note that sink is not on the stack */
  if ( bos != NULL ) {
    for ( i = tos; i != bos; i = i -> bNext ) {
      a = i -> first;
      while ( i -> excess > 0 ) {
	if (( cap[a - G->arcs] == 0 ) && ( a -> resCap > 0 )) {
	  if (a->resCap < i->excess)
	    delta = a->resCap;
	  else
	    delta = i->excess;
	  a -> resCap -= delta;
	  a -> rev -> resCap += delta;
	  i -> excess -= delta;
	  a -> head -> excess += delta;
	}
	a++;
      }
    }
    /* now do the bottom */
    i = bos;
    a = i -> first;
    while ( i -> excess > 0 ) {
      if (( cap[a - G->arcs] == 0 ) && ( a -> resCap > 0 )) {
	if (a->resCap < i->excess)
	  delta = a->resCap;
	else
	  delta = i->excess;
	a -> resCap -= delta;
	a -> rev -> resCap += delta;
	i -> excess -= delta;
	a -> head -> excess += delta;
      }
      a++;
    }
  }
}


/* gap relabeling */

int gap (bucket *emptyB, Graph G)
{

  bucket *l;
  node  *i; 
  long  r;           /* index of the bucket before l  */
  int   cc;          /* cc = 1 if no nodes with positive excess before
		      the gap */

  G->gapCnt ++;
  r = ( emptyB - G->buckets ) - 1;

  /* set labels of nodes beyond the gap to "infinity" */
  for ( l = emptyB + 1; l <= G->buckets + G->dMax; l ++ ) {
    /* this does nothing for high level selection 
    for (i = l -> firstActive; i != sentinelNode; i = i -> bNext) {
      i -> d = n;
      gNodeCnt++;
    }
    l -> firstActive = sentinelNode;
    */

    for ( i = l -> firstInactive; i != G->sentinelNode; i = i -> bNext ) {
      i -> d = n;
      G->gNodeCnt ++;
    }

    l -> firstInactive = G->sentinelNode;
  }

  cc = ( G->aMin > r ) ? 1 : 0;

  G->dMax = r;
  G->aMax = r;

  return ( cc );

}

/*--- relabelling node i */

long relabel (node *i, Graph G)
{

  node  *j;
  long  minD;     /* minimum d of a node reachable from i */
  arc   *minA;    /* an arc which leads to the node with minimal d */
  arc   *a;

  assert(i->excess > 0);

  G->relabelCnt++;
  G->workSinceUpdate += BETA;

  i->d = minD = n;
  minA = NULL;

  /* find the minimum */
  forAllArcs(i,a, G) {
    G->workSinceUpdate++;
    if (a -> resCap > 0) {
      j = a -> head;
      if (j->d < minD) {
	minD = j->d;
	minA = a;
      }
    }
  }

  minD++;
      
  if (minD < n) {

    i->d = minD;
    i->current = minA;

    if (G->dMax < minD) G->dMax = minD;

  } /* end of minD < n */
      
  return ( minD );

} /* end of relabel */


/* discharge: push flow out of i until i becomes inactive */

void discharge (node *i, Graph G)
{

  node  *j;                 /* sucsessor of i */
  long  jD;                 /* d of the next bucket */
  bucket *lj;               /* j's bucket */
  bucket *l;                /* i's bucket */
  arc   *a;                 /* current arc (i,j) */
  cType  delta;
  arc *stopA;

  assert(i->excess > 0);
  assert(i != G->sink);
  do {

    jD = i->d - 1;
    l = G->buckets + i->d;

    /* scanning arcs outgoing from  i  */
    for (a = i->current, stopA = (i+1)->first; a != stopA; a++) {
      if (a -> resCap > 0) {
	j = a -> head;

	if (j->d == jD) {
	  G->pushCnt ++;
	  if (a->resCap < i->excess)
	    delta = a->resCap;
	  else
	    delta = i->excess;
	  a->resCap -= delta;
	  a->rev->resCap += delta;

	  if (j != G->sink) {

	    lj = G->buckets + jD;

	    if (j->excess == 0) {
	      /* remove j from inactive list */
	      iDelete(lj,j, G);
	      /* add j to active list */
	      aAdd(lj,j, G);
	    }
	  }

	  j -> excess += delta;
	  i -> excess -= delta;
	  
	  if (i->excess == 0) break;

	} /* j belongs to the next bucket */
      } /* a  is not saturated */
    } /* end of scanning arcs from  i */

    if (a == stopA) {
      /* i must be relabeled */
      relabel (i, G);

      if (i->d == n) break;
      if ((l -> firstActive == G->sentinelNode) && 
	  (l -> firstInactive == G->sentinelNode)
	  )
	gap (l, G);

      if (i->d == n) break;
    }
    else {
      /* i no longer active */
      i->current = a;
      /* put i on inactive list */
      iAdd(l,i, G);
      break;
    }
  } while (1);
}


// go from higher to lower buckets, push flow
void wave( Graph G) {

  node   *i;
  bucket  *l;

  for (l = G->buckets + G->aMax; l > G->buckets; l--) {
    for (i = l->firstActive; i != G->sentinelNode; i = l->firstActive) {
      aRemove(l,i);

      assert(i->excess > 0);
      discharge (i, G);

    }
  }
}


/* first stage  -- maximum preflow*/

void stageOne ( Graph G )

{

  node   *i;
  bucket  *l;             /* current bucket */


#if defined(INIT_UPDATE) || defined(OLD_INIT) || defined(WAVE_INIT)
  globalUpdate ();
#endif

  G->workSinceUpdate = 0;

#ifdef WAVE_INIT
  wave();
#endif  

  /* main loop */
  while ( G->aMax >= G->aMin ) {
    l = G->buckets + G->aMax;
    i = l->firstActive;

    if (i == G->sentinelNode)
      G->aMax--;
    else {
      aRemove(l,i);

      assert(i->excess > 0);
      discharge (i, G);

      if (G->aMax < G->aMin)
	break;

      /* is it time for global update? */
      if (G->workSinceUpdate * globUpdtFreq > nm) {
	globalUpdate (G);
	G->workSinceUpdate = 0;
      }
      // check wether the sink is still the neighbor of the source 
      // in the tree or if it was changed by another Thread
#ifdef SHORTFLOW
      if (G->sink != G->nodes + tree[Index(G->source, G)])
	return ;
#endif
    }
  } /* end of the main loop */
    
   G->flow = G->sink -> excess;
} 

void check_solution(Graph g)
{
  node *i;
  arc *a;
  excessType sum;
  bucket *l;

  /* check if you have a flow (pseudoflow) */
  /* check arc flows */
  //printf("before tree - t%d step %d\n", id, step);
  forAllNodes(i, g) {
    forAllArcs(i,a, g) {
      if (cap[a - g->arcs] > 0) { /* original arc */        
	if ((a->resCap + a->rev->resCap != cap[a - g->arcs]) 
	    || (a->resCap < 0)
	    || (a->rev->resCap < 0)) {
	  printf("ERROR: bad arc flow\n");
	  exit(2);
	}
      }
    }
  }
  /* check conservation */
  stageTwo(g);
  forAllNodes(i, g)
    if ((i != g->source) && (i != g->sink)) {
      if (i->excess != 0) {
	printf("ERROR: nonzero node excess\n");
	exit(2);
      }
      sum = 0;
      forAllArcs(i,a, g) {
	if (cap[a - g->arcs] > 0) /* original arc */
	  sum -= cap[a - g->arcs] - a->resCap;
	else
	  sum += a->resCap;
      }
      
      if (i->excess != sum) {
	printf("ERROR: conservation constraint violated\n");
	exit(2);
      }
    }
  /* check if mincut is saturated */
  g->aMax = g->dMax = 0;
  for (l = g->buckets; l < g->buckets + n; l++) {
    l->firstActive = g->sentinelNode;
    l->firstInactive = g->sentinelNode;
  }
  globalUpdate(g);
  if (g->source->d < n) {
    printf("ERROR: the solution is not optimal\n");
    exit(2);
  }
  printf("c Solution checks (feasible and optimal)\n");
}

void print_stat(Graph G)
{
      printf ("c pushes:      %10ld\n", G->pushCnt);
      printf ("c relabels:    %10ld\n", G->relabelCnt);
      printf ("c updates:     %10ld\n", G->updateCnt);
      printf ("c gaps:        %10ld\n", G->gapCnt);
      printf ("c gap nodes:   %10ld\n", G->gNodeCnt);
      printf ("c\n");
}

void print_graph(Graph G)
{
  node *i;
  arc *a;
  
  forAllNodes(i, G) {
    long ni = nNode(i, G);
    printf("n=%ld\n", ni);
    forAllArcs(i, a, G) {
      printf("{%ld, %ld}, ", ni, nNode( a -> head, G ));
    }
    printf("\n");
  }   
}  

int RandomInteger (int low, int high);

typedef struct neighbor {
  long node;
  struct neighbor *next;
} * neighbor;  



int main (int argc, char *argv[])
{
  node *j;
  int  cc;

  if (argc > 4 || argc < 3) {
    printf("Usage: %s nthreads inputfile [update frequency]\n", argv[0]);
    exit(1);
  }

  if (argc < 4)
    seed = 1;
  else seed = (unsigned) atoi(argv[3]);

  if (argc != 5)
    globUpdtFreq = GLOB_UPDT_FREQ;
  else
    globUpdtFreq = (float) atof(argv[4]);

  int NTHREADS = (int) atoi(argv[1]);

  // printf("c\nc hi_pr version 3.6\n");
  // printf("c Copyright C by IG Systems, igsys@eclipse.net\nc\n");

  // Allocate and read one copy of the graph for each thread
  Graph *G = (Graph *) calloc(NTHREADS, sizeof(Graph));
  if (G == NULL) {
      fprintf ( stderr, "Allocation error\n"); exit ( 1 );
  }
  int k;
  long max_degree = -1;
  for (k=0; k<NTHREADS; k++) {
    G[k] = (Graph ) malloc(sizeof(struct graph));
    if (G[k] == NULL) {
      fprintf ( stderr, "Allocation error\n"); exit ( 1 );
    }
    parse(argv[2], &n, &m, &(G[k]->nodes), &(G[k]->arcs), &cap, &nMin );
    cc = allocDS(G[k]);
    if ( cc ) { fprintf ( stderr, "Allocation error\n"); exit ( 1 ); }
    // weighted degree computation (capacities of the trivial cuts - jc 03/2011)
    node *v;
    arc *a;
    forAllNodes(v, G[k]) {
      v->degree = 0;
      forAllArcs(v, a, G[k]) {
	v->degree += a->resCap;
      }
      if (v->degree > max_degree)
	max_degree = v->degree;
    }
  }
  // end of allocate and read graphs 

  // initial permutation of the nodes
  long *permut = (long *)calloc(n, sizeof(long));
  if (permut == NULL)
    {
      fprintf ( stderr, "Allocation error\n"); exit ( 1 );
    }
  int i;
  for (i = 0; i < n; i++)
    permut[i] = i;
#ifndef DONTPERMUTE
  srand (seed); 
  for (i = 0; i < n-1; i++)
    {
      long s_rand = RandomInteger(i, n-1);
      long tmp = permut[s_rand];
      permut[s_rand] = permut[i];
      permut[i] = tmp;
    } 
#endif

#ifdef DEGREE_HEURISTIC
  long *count = (long *)calloc(max_degree+1, sizeof(long));
  long *permut2 = (long *)calloc(n, sizeof(long));
  if (count == NULL || permut2 == NULL)
    {
      fprintf ( stderr, "Allocation error\n"); exit ( 1 );
    }
  // counting sort the nodes by degree
  for (i=0; i<=max_degree; i++) count[i] = 0;
  for (i=0; i<n; i++) count[ (G[0]->nodes + permut[i])->degree]++;
  for (i=1; i<=max_degree; i++) count[i] += count[i-1];
  for (i=0; i<n; i++) 
#ifdef DECREASING
    permut2[n - count[(G[0]->nodes+permut[i])->degree]--] = permut[i];
#else
    permut2[--count[(G[0]->nodes+permut[i])->degree]] = permut[i];
#endif

#ifdef VERBOSE
  for (i=0; i<n; i++)
    printf("c Node: %ld Degree: %ld\n", permut2[i]+1, (G[0]->nodes+permut2[i])->degree); 
#endif
  free(permut);
  permut = permut2;
  free(count);
#endif

  // optimization with neighbor lists - jaime 29/10/2010
  neighbor st_neighbors = (neighbor) calloc(n+1, sizeof(struct neighbor));
  neighbor *neighbors = (neighbor*) calloc(n, sizeof(neighbor));
  if (st_neighbors == NULL || neighbors == NULL)
    {
      printf("cannot obtain enough memory.\n");
      exit(1);
    }

  printf("c nodes:       %10ld\nc arcs:        %10ld\nc\n", n, m);

  timer1();
  double topenstart = omp_get_wtime();

  // Tree initialization:
  // the tree is represented by a vector tree[] such that each node points to its parent
  // plus linked lists where each node points to its neighbors that have not yet been 
  // separated by a cut.
  tree = (long *) calloc(n, sizeof(long)); 
  tree_weights = (double *) calloc (n, sizeof(double)); 
  if (tree == NULL || tree_weights == NULL) 
    {
      printf("can't obtain enough memory to solve this problem.\n");
      exit(1);
    }
  tree[permut[0]] = -1; 
  tree_weights[permut[0]] = -1.0;
  st_neighbors[0].node = permut[0];
  neighbors[permut[0]] = &st_neighbors[1];
  for (i=1; i < n; ++i) 
    {
      tree[permut[i]] = permut[0];
      tree_weights[permut[i]] = -1.0;
      st_neighbors[i].node = permut[i];
      st_neighbors[i].next = &st_neighbors[i+1];      
      neighbors[permut[i]] = NULL;
    }
  st_neighbors[n-1].next = NULL;

  long success=0, unsuccess=0;

  //--------------------------------------------------------------------------
  //   MAIN LOOP
  //---------------------------------------------------------------------------
  // gusfield iteration starts here 
  // , i, a, l, cap, sum
#pragma omp parallel for default(none) \
  num_threads(NTHREADS) \
  private(step) shared(j, tree,tree_weights,G,n, success, unsuccess, \
		       permut, neighbors) 
  // schedule(guided) 
  // reduction(+:success, unsuccess)	//lar
  for (step=1; step < n; step++)
    {
      int id = omp_get_thread_num();
      int done = 0;

      // set the source
      long  source = permut[step];
      G[id]->source = G[id]->nodes + source; 

     do { // do until the cut is usefull for this source 

       // set the sink 
       G[id]->sink = G[id]->nodes +  tree[source]; 
      
       // run and time MaxFlow (stageOne)
       init( G[id] );
       stageOne ( G[id] );

#ifdef VERBOSE
       printf("c Thread %d: source=%ld sink=%ld flow=%.1lf\n", id, source+1, tree[source]+1, G[id]->flow);
#endif 

       // adjust labels
       globalUpdate(G[id]); 

       // verify optimality
       if (G[id]->source->d < n) {
	 printf("ERROR: the solution is not optimal\n");
	 exit(2);
       }

#ifdef CHECK_SOLUTION
      check_solution(G[id]);
#endif
#ifdef PRINT_STAT
      print_stat(G[id]);
#endif

      // adjust the tree - CRITICAL REGION
#pragma omp critical //(lar)
{ 
      long target = tree[source];
      if (G[id]->sink == G[id]->nodes + target) {
	done = 1;
	success++;
	tree_weights[source] = G[id]->flow; 
	
	// printf("c >>> source %ld  target %ld\n", source+1,target+1);

#ifndef NOTRIVIALCUT
	if (G[id]->source->degree > G[id]->flow)
	  {
#endif
	// for all target's neighbor
	    neighbor nb, prev_nb, next_nb;
	    for (nb = neighbors[target], prev_nb = NULL; nb != NULL; nb = next_nb)
	      {
		next_nb = nb->next;
		j = G[id]->nodes + nb->node; 
		if (j->d >= n)  // j is a node on the source side
		  {
		    // remove nb from target's list
		    if (prev_nb == NULL)
		      neighbors[target] = nb->next;
		    else
		      prev_nb->next = nb->next;
		    // move nb from target list to source list
		    if (nb->node != source)
		      {
			nb->next = neighbors[source];
			neighbors[source] = nb;
			tree[nb->node] = source;   // nb on the 'source' side
		      }
		    // in both cases above, prev_nb remains unchanged
		    // because nb was removed from the list
		  } else prev_nb = nb;
	      }
#ifndef NOTRIVIALCUT 
	  } // end if d(v) > flow
#endif
      } else { // cannot use the cut and it will compute another one
	unsuccess++;
#ifdef VERBOSE
	printf("c fail thread=%d source=%ld sink=%d new_sink=%ld\n",id, source+1, Index(G[id]->sink, G[id])+1, tree[source]+1);
#endif
      }
 } // end of omp critical 
     } while(!done);
    } // end of gus main loop (for)
  // -- End of parallel for ------------------------------------------------
  
  // begin: print the flow equivalent tree 
  printf("c flow equivalent tree\n");
  printf("p cut %ld %ld\n", n, n-1);
  for (step = 1; step < n; ++step)
    printf("a %ld %ld %lf\n", permut[step]+1, tree[permut[step]]+1, tree_weights[permut[step]]);
  // end: print tree

  double topenfinish = omp_get_wtime();
  timer2();
  printf("c timer: real: %.7lf user: %.7lf sys: %.7lf\n", topenfinish-topenstart, getUserTime(), getSystemTime()); 
  printf("c success: %ld unsuccess: %ld\n",success, unsuccess);

  exit(0);
}

int RandomInteger (int low, int high)
{
    int k;
    double d;
    d = (double) rand () / ((double) RAND_MAX + 1);
    k = d * (high - low + 1);
    return low + k;
}
