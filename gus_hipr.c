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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

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
node   *nodes;               /* array of nodes */
arc    *arcs;                /* array of arcs */
bucket *buckets;             /* array of buckets */
cType  *cap;                 /* array of capacities */
node   *source;              /* source node pointer */
node   *sink;                /* sink node pointer */
//node   **queue;              /* queue for BFS */
//node   **qHead, **qTail, **qLast;     /* queue pointers */
long   dMax;                 /* maximum label */
long   aMax;                 /* maximum actie node label */
long   aMin;                 /* minimum active node label */
double flow;                 /* flow value */
long pushCnt  = 0;           /* number of pushes */
long relabelCnt   = 0;       /* number of relabels */
long updateCnt    = 0;       /* number of updates */
long gapCnt   = 0;           /* number of gaps */
long gNodeCnt = 0;           /* number of nodes after gap */  
float t, t2;                 /* for saving times */
node   *sentinelNode;        /* end of the node list marker */
arc *stopA;                  /* used in forAllArcs */
long workSinceUpdate=0;      /* the number of arc scans since last update */
float globUpdtFreq;          /* global update frequency */

// begin: variables declarations (jc)
long *tree; 
double *tree_weights;
long step; 
// end: variables (jc) 
int RandomInteger (int low, int high);

/* macros */

#define forAllNodes(i) for ( i = nodes; i != sentinelNode; i++ )
#define forAllArcs(i,a) for (a = i->first, stopA = (i+1)->first; a != stopA; a++)

#define nNode( i ) ( (i) - nodes + nMin )
// Index(): from pointer to index (jc)
#define Index( i ) ( (i) - nodes )
#define nArc( a )  ( ( a == NULL )? -1 : (a) - arcs )

#define min( a, b ) ( ( (a) < (b) ) ? a : b )

long i_dist;

#define aAdd(l,i)\
{\
  i->bNext = l->firstActive;\
  l->firstActive = i;\
  i_dist = i->d;\
  if (i_dist < aMin)\
    aMin = i_dist;\
  if (i_dist > aMax)\
    aMax = i_dist;\
  if (dMax < aMax)\
    dMax = aMax;\
}

/* i must be the first element */
#define aRemove(l,i)\
{\
  l->firstActive = i->bNext;\
}

node *i_next, *i_prev;
#define iAdd(l,i)\
{\
  i_next = l->firstInactive;\
  i->bNext = i_next;\
  i->bPrev = sentinelNode;\
  i_next->bPrev = i;\
  l->firstInactive = i;\
}

#define iDelete(l,i)\
{\
  i_next = i->bNext;\
  if (l->firstInactive == i) {\
    l->firstInactive = i_next;\
    i_next->bPrev = sentinelNode;\
  }\
  else {\
    i_prev = i->bPrev;\
    i_prev->bNext = i_next;\
    i_next->bPrev = i_prev;\
  }\
}

/* allocate datastructures, initialize related variables */

int allocDS( )

{

  nm = ALPHA * n + m;
  /*
  queue = (node**) calloc ( n, sizeof (node*) );
  if ( queue == NULL ) return ( 1 );
  qLast = queue + n - 1;
  qInit();
  */
  buckets = (bucket*) calloc ( n+2, sizeof (bucket) );
  if ( buckets == NULL ) return ( 1 );

  sentinelNode = nodes + n;
  sentinelNode->first = arcs + 2*m;

  return ( 0 );

} /* end of allocate */


void init( )

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

  // initialize excesses

  forAllNodes(i) {
    i->excess = 0;
    i->current = i->first;
    forAllArcs(i, a)
      {
	a->resCap = cap[a-arcs];
	// a->rev->resCap = cap[a-arcs]; // for undirected graphs, did not work. jc
      }
  }

  for (l = buckets; l <= buckets + n-1; l++) {
    l -> firstActive   = sentinelNode;
    l -> firstInactive  = sentinelNode;
  }
    
  overflowDetected = 0;
#ifdef EXCESS_TYPE_LONG
  testExcess = 0;
  forAllArcs(source,a) {
    if (a->head != source) {
      testExcess += a->resCap;
    }
  }
  if (testExcess > MAXLONG) {
    printf("c WARNING: excess overflow. See README for details.\nc\n");
    overflowDetected = 1;
  }
#endif
#ifdef OLD_INIT
  source -> excess = MAXLONG;
#else
  if (overflowDetected) {
    source -> excess = MAXLONG;
  }
  else {
    source->excess = 0;
    forAllArcs(source,a) {
      if (a->head != source) {
	pushCnt ++;
	delta = a -> resCap;
	a -> resCap -= delta;
	(a -> rev) -> resCap += delta;
	a->head->excess += delta;
      }
    }
  }

  /*  setup labels and buckets */
  l = buckets + 1;
    
  aMax = 0;
  aMin = n;
    
  forAllNodes(i) {
    if (i == sink) {
      i->d = 0;
      iAdd(buckets,i);
      continue;
    }
    if ((i == source) && (!overflowDetected)) {
      i->d = n;
    }
    else
      i->d = 1;
    if (i->excess > 0) {
      /* put into active list */
      aAdd(l,i);
    }
    else { /* i -> excess == 0 */
      /* put into inactive list */
      if (i->d < n)
	iAdd(l,i);
    }
  }
  dMax = 1;
#endif

  //  dMax = n-1;
  //  flow = 0.0;

} /* end of init */

void checkMax()

{
  bucket *l;

  for (l = buckets + dMax + 1; l < buckets + n; l++) {
    assert(l->firstActive == sentinelNode);
    assert(l->firstInactive == sentinelNode);
  }
}

/* global update via backward breadth first search from the sink */

void globalUpdate ()

{

  node  *i, *j;       /* node pointers */
  arc   *a;           /* current arc pointers  */
  bucket *l, *jL;          /* bucket */
  long curDist, jD;
  long state;


  updateCnt ++;

  /* initialization */

  forAllNodes(i)
    i -> d = n;
  sink -> d = 0;

  dMax = n; // (jc - essa linha nao existia - resolveu o erro do dblcyc.1024 )
  for (l = buckets; l <= buckets + dMax; l++) {
    l -> firstActive   = sentinelNode;
    l -> firstInactive  = sentinelNode;
  }

  dMax = aMax = 0;
  aMin = n;

  /* breadth first search */

  // add sink to bucket zero

  iAdd(buckets, sink);
  for (curDist = 0; 1; curDist++) {

    state = 0;
    l = buckets + curDist;
    jD = curDist + 1;
    jL = l + 1;
    /*
    jL -> firstActive   = sentinelNode;
    jL -> firstInactive  = sentinelNode;
    */

    if ((l->firstActive == sentinelNode) && 
	(l->firstInactive == sentinelNode))
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

       if (i == sentinelNode) {
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
      forAllArcs(i,a) {
	if (a->rev->resCap > 0 ) {
	  j = a->head;
	  if (j->d == n) {
	    j->d = jD;
	    j->current = j->first;
	    if (jD > dMax) dMax = jD;
	    
	    if (j->excess > 0) {
	      /* put into active list */
	      aAdd(jL,j);
	    }
	    else {
	      /* put into inactive list */
	      iAdd(jL,j);
	    }
	  }
	}
      } /* node i is scanned */ 
    }
  }

} /* end of global update */

/* second stage -- preflow to flow */
void stageTwo ( )
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
  forAllNodes(i) {
    forAllArcs(i,a)
      if ( a -> head == i ) {
	a -> resCap = cap[a - arcs];
      }
  }

  /* initialize */
  tos = bos = NULL;
  forAllNodes(i) {
    i -> d = WHITE;
    //    buckets[i-nodes].firstActive = NULL;
    buckets[i-nodes].firstActive = sentinelNode;
    i -> current = i -> first;
  }

  /* eliminate flow cycles, topologicaly order vertices */
  forAllNodes(i)
    if (( i -> d == WHITE ) && ( i -> excess > 0 ) &&
	( i != source ) && ( i != sink )) {
      r = i;
      r -> d = GREY;
      do {
	for ( ; i->current != (i+1)->first; i->current++) {
	  a = i -> current;
	  if (( cap[a - arcs] == 0 ) && ( a -> resCap > 0 )) { 
	    j = a -> head;
	    if ( j -> d == WHITE ) {
	      /* start scanning j */
	      j -> d = GREY;
	      buckets[j-nodes].firstActive = i;
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
	  if ( i != source ) {
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
	    i = buckets[i-nodes].firstActive;
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
	if (( cap[a - arcs] == 0 ) && ( a -> resCap > 0 )) {
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
      if (( cap[a - arcs] == 0 ) && ( a -> resCap > 0 )) {
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

int gap (emptyB)
     bucket *emptyB;

{

  bucket *l;
  node  *i; 
  long  r;           /* index of the bucket before l  */
  int   cc;          /* cc = 1 if no nodes with positive excess before
		      the gap */

  gapCnt ++;
  r = ( emptyB - buckets ) - 1;

  /* set labels of nodes beyond the gap to "infinity" */
  for ( l = emptyB + 1; l <= buckets + dMax; l ++ ) {
    /* this does nothing for high level selection 
    for (i = l -> firstActive; i != sentinelNode; i = i -> bNext) {
      i -> d = n;
      gNodeCnt++;
    }
    l -> firstActive = sentinelNode;
    */

    for ( i = l -> firstInactive; i != sentinelNode; i = i -> bNext ) {
      i -> d = n;
      gNodeCnt ++;
    }

    l -> firstInactive = sentinelNode;
  }

  cc = ( aMin > r ) ? 1 : 0;

  dMax = r;
  aMax = r;

  return ( cc );

}

/*--- relabelling node i */

long relabel (i)

node *i;   /* node to relabel */

{

  node  *j;
  long  minD;     /* minimum d of a node reachable from i */
  arc   *minA;    /* an arc which leads to the node with minimal d */
  arc   *a;

  assert(i->excess > 0);

  relabelCnt++;
  workSinceUpdate += BETA;

  i->d = minD = n;
  minA = NULL;

  /* find the minimum */
  forAllArcs(i,a) {
    workSinceUpdate++;
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

    if (dMax < minD) dMax = minD;

  } /* end of minD < n */
      
  return ( minD );

} /* end of relabel */


/* discharge: push flow out of i until i becomes inactive */

void discharge (i)

node  *i;

{

  node  *j;                 /* sucsessor of i */
  long  jD;                 /* d of the next bucket */
  bucket *lj;               /* j's bucket */
  bucket *l;                /* i's bucket */
  arc   *a;                 /* current arc (i,j) */
  cType  delta;
  arc *stopA;

  assert(i->excess > 0);
  assert(i != sink);
  do {

    jD = i->d - 1;
    l = buckets + i->d;

    /* scanning arcs outgoing from  i  */
    for (a = i->current, stopA = (i+1)->first; a != stopA; a++) {
      if (a -> resCap > 0) {
	j = a -> head;

	if (j->d == jD) {
	  pushCnt ++;
	  if (a->resCap < i->excess)
	    delta = a->resCap;
	  else
	    delta = i->excess;
	  a->resCap -= delta;
	  a->rev->resCap += delta;

	  if (j != sink) {

	    lj = buckets + jD;

	    if (j->excess == 0) {
	      /* remove j from inactive list */
	      iDelete(lj,j);
	      /* add j to active list */
	      aAdd(lj,j);
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
      relabel (i);

      if (i->d == n) break;
      if ((l -> firstActive == sentinelNode) && 
	  (l -> firstInactive == sentinelNode)
	  )
	gap (l);

      if (i->d == n) break;
    }
    else {
      /* i no longer active */
      i->current = a;
      /* put i on inactive list */
      iAdd(l,i);
      break;
    }
  } while (1);
}


// go from higher to lower buckets, push flow
void wave() {

  node   *i;
  bucket  *l;

  for (l = buckets + aMax; l > buckets; l--) {
    for (i = l->firstActive; i != sentinelNode; i = l->firstActive) {
      aRemove(l,i);

      assert(i->excess > 0);
      discharge (i);

    }
  }
}


/* first stage  -- maximum preflow*/

void stageOne ( )

{

  node   *i;
  bucket  *l;             /* current bucket */


#if defined(INIT_UPDATE) || defined(OLD_INIT) || defined(WAVE_INIT)
  globalUpdate ();
#endif

  workSinceUpdate = 0;

#ifdef WAVE_INIT
  wave();
#endif  

  /* main loop */
  while ( aMax >= aMin ) {
    l = buckets + aMax;
    i = l->firstActive;

    if (i == sentinelNode)
      aMax--;
    else {
      aRemove(l,i);

      assert(i->excess > 0);
      discharge (i);

      if (aMax < aMin)
	break;

      /* is it time for global update? */
      if (workSinceUpdate * globUpdtFreq > nm) {
	globalUpdate ();
	workSinceUpdate = 0;
      }

    }
    
  } /* end of the main loop */
    
  flow = sink -> excess;

} 


int main (argc, argv)

     int   argc;
     char *argv[];

{
#if (defined(PRINT_FLOW) || defined(CHECK_SOLUTION))
  node *i;
  arc *a;
#endif

#ifdef PRINT_FLOW
  long ni, na;
#endif

  node *j;
  int  cc;
#ifdef CHECK_SOLUTION
  excessType sum;
  bucket *l;
#endif


  if (argc > 3 || argc < 2) {
    printf("Usage: %s inputfile [update frequency]\n", argv[0]);
    exit(1);
  }

  if (argc != 3)
    globUpdtFreq = GLOB_UPDT_FREQ;
  else
    globUpdtFreq = (float) atof(argv[2]);

  // printf("c\nc hi_pr version 3.6\n");
  // printf("c Copyright C by IG Systems, igsys@eclipse.net\nc\n");

  parse(argv[1], &n, &m, &nodes, &arcs, &cap, &nMin );

  // printf("c nodes:       %10ld\nc arcs:        %10ld\nc\n", n, m);

  cc = allocDS();
  if ( cc ) { fprintf ( stderr, "Allocation error\n"); exit ( 1 ); }

  timer1();

  // initialization for gusfield (jc)
  tree = (long *)calloc(n, sizeof(long)); // jc
  tree_weights = (double *)calloc (n, sizeof(double)); // jc
  if (tree == NULL || tree_weights == NULL) // jc
    {
      printf("can't obtain enough memory to solve this problem.\n");
      exit(1);
    }

  // initial permutation of the nodes
  long *permut = (long *)calloc(n, sizeof(long));
  int i;
  for (i = 0; i < n; i++)
    permut[i] = i;
#ifndef DONTPERMUTE
  srand (time (NULL)); 
  for (i = 0; i < n-1; i++)
    {
      long s_rand = RandomInteger(i, n-1);
      long tmp = permut[s_rand];
      permut[s_rand] = permut[i];
      permut[i] = tmp;
    } 
#endif

  tree[permut[0]] = -1; 
  tree_weights[permut[0]] = -1.0;
  for (i=1; i < n; ++i) 
    {
      tree[permut[i]] = permut[0];
      tree_weights[permut[i]] = -1.0;
    }

  // gusfield iteration starts here
  for (step=1; step < n; step++)
    {
      //printf("source = %ld   sink = %ld\n", step+1, tree[step]+1);
      long  i_source = permut[step];
      source = nodes + i_source; 
      sink = nodes +  tree[i_source]; 
      t = timer();
      t2 = t;
      init();
      stageOne ( );
      t2 = timer() - t2;
      // printf ("c flow:       %12.01f\n", flow);
#ifndef CUT_ONLY
      stageTwo ( );
      t = timer() - t;
      // printf ("c time:        %10.2f\n", t);
#endif
      // printf ("c cut tm:      %10.2f\n", t2);
#ifdef CHECK_SOLUTION
      /* check if you have a flow (pseudoflow) */
      /* check arc flows */
      forAllNodes(i) {
	forAllArcs(i,a) {
	  if (cap[a - arcs] > 0) /* original arc */
	    if ((a->resCap + a->rev->resCap != cap[a - arcs]) 
		|| (a->resCap < 0)
		|| (a->rev->resCap < 0)) {
	      printf("ERROR: bad arc flow\n");
	      exit(2);
	    }
	}
      }
      /* check conservation */
      forAllNodes(i)
	if ((i != source) && (i != sink)) {
#ifdef CUT_ONLY
	  if (i->excess < 0) {
	    printf("ERROR: nonzero node excess\n");
	    exit(2);
	  }
#else
	  if (i->excess != 0) {
	    printf("ERROR: nonzero node excess\n");
	    exit(2);
	  }
#endif
	  sum = 0;
	  forAllArcs(i,a) {
	    if (cap[a - arcs] > 0) /* original arc */
	      sum -= cap[a - arcs] - a->resCap;
	    else
	      sum += a->resCap;
	  }

	  if (i->excess != sum) {
	    printf("ERROR: conservation constraint violated\n");
	    exit(2);
	  }
	}
      /* check if mincut is saturated */
      aMax = dMax = 0;
      for (l = buckets; l < buckets + n; l++) {
	l->firstActive = sentinelNode;
	l->firstInactive = sentinelNode;
      }
      globalUpdate();
      if (source->d < n) {
	printf("ERROR: the solution is not optimal\n");
	exit(2);
      }
      printf("c\nc Solution checks (feasible and optimal)\nc\n");
#endif

#ifdef PRINT_STAT
      printf ("c pushes:      %10ld\n", pushCnt);
      printf ("c relabels:    %10ld\n", relabelCnt);
      printf ("c updates:     %10ld\n", updateCnt);
      printf ("c gaps:        %10ld\n", gapCnt);
      printf ("c gap nodes:   %10ld\n", gNodeCnt);
      printf ("c\n");
#endif

#ifdef PRINT_FLOW
      printf ("c flow values\n");
      forAllNodes(i) {
	ni = nNode(i);
	forAllArcs(i,a) {
	  na = nArc(a);
	  if ( cap[na] > 0 )
	    printf ( "f %7ld %7ld %12ld\n",
		     ni, nNode( a -> head ), cap[na] - ( a -> resCap )
		     );
	}
      }
      printf("c\n");
#endif

      globalUpdate(); 
      /*
      printf ("c nodes on the sink side\n");
      forAllNodes(j)
	if (j->d < n)
	  printf("c %ld ", nNode(j));
      printf("\n");
      */

      // adjust the tree
      tree_weights[i_source] = flow; 
      forAllNodes(j) 
	if (j->d >= n) // j is a node on the source side
	  if (tree_weights[Index(j)] < 0.0 && (tree[Index(j)] == tree[i_source]))
	    tree[Index(j)] = i_source;   // now nNode(j) is adjacent to 'step'
    } // gus iteration ends here (jc)

  timer2();

  // begin: print the flow equivalent tree (jc)
  printf("c flow equivalent tree\n");
  printf("p cut %ld %ld\n", n, n-1);
  for (step = 1; step < n; ++step)
    printf("a %ld %ld %.1lf\n", permut[step]+1, tree[permut[step]]+1, tree_weights[permut[step]]);
  printf("c timer: real: %.7lf user: %.7lf sys: %.7lf\n", getRealTime(), getUserTime(), getSystemTime()); 
  // end: print tree

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
