#include <unistd.h>
#include <sys/wait.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <limits.h>
#include "misc.h"
#include "tinygraph.h"
#include "graph.h"
#include "heap.h"
#include "blant.h"
#include "queue.h"
#include "multisets.h"
#include "sorts.h"
#include "blant-window.h"
#include "blant-output.h"
#include "blant-utils.h"
#include "rand48.h"
Boolean _earlyAbort; // Can be set true by anybody anywhere, and they're responsible for producing a warning as to why

static int *_pairs, _numNodes, _numEdges, _maxEdges=1024, _seed = -1; // -1 means "not initialized"
char **_nodeNames, _supportNodeNames = true;
Boolean _child; // are we a child process?

char * _sampleFileName;

// _k is the global variable storing k; _Bk=actual number of entries in the canon_map for given k.
unsigned int _k;
unsigned int _Bk, _k_small;

int _alphaList[MAX_CANONICALS];
int _numCanon, _numSamples;
Gint_type _canonList[MAX_CANONICALS]; // map ordinals to integer representation of the canonical
SET *_connectedCanonicals; // the SET of canonicals that are connected.
int _numConnectedCanon;
int _numConnectedComponents;
int *_componentSize;

int _numOrbits, _orbitList[MAX_CANONICALS][MAX_K]; // map from [ordinal][canonicalNode] to orbit ID.
int _orbitCanonMapping[MAX_ORBITS]; // Maps orbits to canonical (including disconnected)
int _orbitCanonNodeMapping[MAX_ORBITS]; // Maps orbits to canonical (including disconnected)
int *_whichComponent;

// char* _BLANT_DIR;

enum OutputMode _outputMode = undef;
unsigned long int _graphletCount[MAX_CANONICALS];
int **_graphletDistributionTable;
double _g_overcount, _graphletConcentration[MAX_CANONICALS];

enum CanonicalDisplayMode _displayMode = undefined;
enum FrequencyDisplayMode _freqDisplayMode = freq_display_mode_undef;

int _outputMapping[MAX_CANONICALS];

int _orca_orbit_mapping[58];
int _connectedOrbits[MAX_ORBITS];
int _numConnectedOrbits;

// A bit counter-intuitive: we need to allocate this many vectors each of length [_numNodes],
// and then the degree for node v, graphlet/orbit g is _degreeVector[g][v], NOT [v][g].
// We do this simply because we know the length of MAX_CANONICALS so we pre-know the length of
// the first dimension, otherwise we'd need to get more funky with the pointer allocation.
// Only one of these actually get allocated, depending upon outputMode.
unsigned long int *_graphletDegreeVector[MAX_CANONICALS];
unsigned long int    *_orbitDegreeVector[MAX_ORBITS];
double *_doubleOrbitDegreeVector[MAX_ORBITS];

double *_cumulativeProb;

// number of parallel threads required, and the maximum allowed at one time.
int _JOBS, _MAX_THREADS;

// Here's the actual mapping from non-canonical to canonical, same argument as above wasting memory, and also mmap'd.
// So here we are allocating 256MB x sizeof(short int) = 512MB.
// Grand total statically allocated memory is exactly 1.25GB.
//static short int _K[maxBk] __attribute__ ((aligned (8192)));
short *_K = NULL; // Allocating memory dynamically

/* AND NOW THE CODE */

// return how many nodes found. If you call it with startingNode == 0 then we automatically clear the visited array
static TSET _visited;
static int NumReachableNodes(TINY_GRAPH *g, int startingNode)
{
    if(startingNode == 0) TSetEmpty(_visited);
    TSetAdd(_visited,startingNode);
    unsigned int j, Varray[MAX_K], numVisited = 0;
    int numNeighbors = TSetToArray(Varray, g->A[startingNode]);
    assert(numNeighbors == g->degree[startingNode]);
    for(j=0; j<numNeighbors; j++)if(!TSetIn(_visited,Varray[j])) numVisited += NumReachableNodes(g,Varray[j]);
    return 1+numVisited;
}

static int **_componentList; // list of lists of components, largest to smallest.
static double _totalCombinations, *_combinations, *_probOfComponent;
SET **_componentSet;

void SetBlantDir() {
    char* temp = getenv("BLANT_DIR");
    if (temp != NULL)
	_BLANT_DIR = strdup(temp); // can't assume the string returned by getetv never changes, so copy it.
}

static int InitializeConnectedComponents(GRAPH *G)
{
    static unsigned int v, *Varray, j, i;
    assert(!Varray); // we only can be called once.
    assert(_numConnectedComponents == 0);
    SET *visited = SetAlloc(G->n);
    Varray = Calloc(G->n, sizeof(int));
    _whichComponent = Calloc(G->n, sizeof(int));
    _componentSize = Calloc(G->n, sizeof(int)); // probably bigger than it needs to be but...
    _componentList = Calloc(G->n, sizeof(int*)); // probably bigger...
    _combinations = Calloc(G->n, sizeof(double*)); // probably bigger...
    _probOfComponent = Calloc(G->n, sizeof(double*)); // probably bigger...
    _cumulativeProb = Calloc(G->n, sizeof(double*)); // probably bigger...
    _componentSet = Calloc(G->n, sizeof(SET*));

    int nextStart = 0;
    _componentList[0] = Varray;
    for(v=0; v < G->n; v++) if(!SetIn(visited, v))
    {
	_componentSet[_numConnectedComponents] = SetAlloc(G->n);
	_componentList[_numConnectedComponents] = Varray + nextStart;
	GraphVisitCC(G, v, _componentSet[_numConnectedComponents], Varray + nextStart, _componentSize + _numConnectedComponents);
	SetUnion(visited, visited, _componentSet[_numConnectedComponents]);
	for(j=0; j < _componentSize[_numConnectedComponents]; j++)
	{
	    assert(_whichComponent[Varray[nextStart + j]] == 0);
	    _whichComponent[Varray[nextStart + j]] = _numConnectedComponents;
	}
	nextStart += _componentSize[_numConnectedComponents];
	++_numConnectedComponents;
    }
    assert(nextStart == G->n);

    _totalCombinations = 0.0;
    for(i=0; i< _numConnectedComponents; i++)
    {
	//find the biggest one
	int biggest = i;
	for(j=i+1; j<_numConnectedComponents;j++)
	    if(_componentSize[j] > _componentSize[biggest])
		biggest = j;
	// Now swap the biggest one into position i;
	for(j=0; j < _componentSize[biggest]; j++)
	    _whichComponent[_componentList[biggest][j]] = i;
	int itmp, *pitmp;
	SET * stmp;
	itmp = _componentSize[i];
	_componentSize[i] = _componentSize[biggest];
	_componentSize[biggest] = itmp;
	pitmp = _componentList[i];
	_componentList[i] = _componentList[biggest];
	_componentList[biggest] = pitmp;
    stmp = _componentSet[i];
    _componentSet[i] = _componentSet[biggest];
    _componentSet[biggest] = stmp;
	_combinations[i] = CombinChooseDouble(_componentSize[i], _k);
	_totalCombinations += _combinations[i];
    }

    double cumulativeProb = 0.0;
    for(i=0; i< _numConnectedComponents; i++)
    {
	_probOfComponent[i] =  _combinations[i] / _totalCombinations;
	_cumulativeProb[i] = cumulativeProb + _probOfComponent[i];
	cumulativeProb = _cumulativeProb[i];
	if(cumulativeProb > 1)
	{
	    assert(cumulativeProb - 1.0 < _numConnectedComponents * 1e-15);  // allow some roundoff error
	    cumulativeProb = _cumulativeProb[i] = 1.0;
	}
	//printf("Component %d has %d nodes and probability %lf, cumulative prob %lf\n", i, _componentSize[i], _probOfComponent[i], _cumulativeProb[i]);
    }
    SetFree(visited);
    return _numConnectedComponents;
}

int alphaListPopulate(char *BUF, int *alpha_list, int k) {
	sprintf(BUF, "%s/%s/alpha_list_mcmc%d.txt", _BLANT_DIR, CANON_DIR, k);
    FILE *fp_ord=fopen(BUF, "r");
    if(!fp_ord) Fatal("cannot find %s\n", BUF);
    int numAlphas, i;
    assert(1==fscanf(fp_ord, "%d",&numAlphas));
	assert(numAlphas == _numCanon);
    for(i=0; i<numAlphas; i++) assert(1==fscanf(fp_ord, "%d", &alpha_list[i]));
    fclose(fp_ord);
    return numAlphas;
}

// Compute the degree of the state in the state graph (see Lu&Bressen)
// Given the big graph G, and a set of nodes S (|S|==k), compute the
// degree of the *state* represented by these nodes, which is:
//     Degree(S) = k*|NN2| + (k-1)*|NN1|
// where NN1 is the set of nodes one step outside S that have only 1 connection back into S,
// and   NN2 is the set of nodes one step outside S that have more than 1 connection back into S.
static int StateDegree(GRAPH *G, SET *S)
{
#if PARANOID_ASSERTS
    assert(SetCardinality(S) == _k);
#endif
    int Varray[_k]; // array of elements in S
    SetToArray(Varray, S);

    SET *outSet=SetAlloc(G->n); // the set of nodes in G that are one step outside S, from *anybody* in S
    int connectCount[G->n]; // for each node in the outset, count the number of times it's connected into S
    memset(connectCount, 0, G->n * sizeof(int)); // set them all to zero

    int i, j, numOut = 0;
    for(i=0;i<_k;i++) for(j=0; j<G->degree[Varray[i]]; j++)
    {
	int neighbor = G->neighbor[Varray[i]][j];
	if(!SetIn(S, neighbor))
	{
	    SetAdd(outSet, neighbor);
	    if(connectCount[neighbor] == 0) ++numOut;
	    ++connectCount[neighbor];
	}
    }
#if PARANOID_ASSERTS
    assert(SetCardinality(outSet) == numOut);
#endif
    int outArray[numOut], Sdeg = 0;
    i = SetToArray(outArray, outSet);
    assert(i == numOut);
    for(i=0; i < numOut; i++)
    {
	switch(connectCount[outArray[i]])
	{
	case 0: break;
	case 1: Sdeg += _k-1; break;
	default:  Sdeg += _k; break;
	}
    }
    SetFree(outSet);
    return Sdeg;
}


void getDoubleDegreeArr(double *double_degree_arr, GRAPH *G) {
    int i;

    for (i = 0; i < G->n; ++i) {
        double_degree_arr[i] = (double)G->degree[i];
    }
}

Boolean ProcessGraphlet(GRAPH *G, SET *V, unsigned Varray[], const int k, TINY_GRAPH *g)
{
    Boolean processed = true;
    TinyGraphInducedFromGraph(g, G, Varray);
    Gint_type Gint = TinyGraph2Int(g,k), GintOrdinal=L_K(Gint), j;

#if PARANOID_ASSERTS
    assert(0 <= GintOrdinal && GintOrdinal < _numCanon);
#endif

	if(NodeSetSeenRecently(G, Varray,k) || !SetIn(_windowRep_allowed_ambig_set, GintOrdinal)) processed=false;
	else puts(PrintIndexEntry(Gint, GintOrdinal, Varray, g, k));

    return processed;
}

void SampleGraphletIndexAndPrint(GRAPH* G, int* prev_nodes_array, int prev_nodes_count, double *heur_arr) {
    // i, j, and neigh are just used in for loops in this function
    int i, j, neigh;
    // the tiny_graph is not used in this function. it is only used as a temporary data object as part of ProcessGraphlet (see below)
    TINY_GRAPH *g = TinyGraphAlloc(_k);

    // base case for the recursion: a k-graphlet is found, print it and return
    if (prev_nodes_count == _k) {
        // ProcessGraphlet will create the k-node induced graphlet from prev_nodes_array, and then determine if said graphlet is of a low enough multiplicity (<= multiplicity)
        // ProcessGraphlet will also check that the k nodes you passed it haven't already been printed (although, this system does not work 100% perfectly)
        // ProcessGraphlet will also print the nodes as output if the graphlet passes all checks
        ProcessGraphlet(G, NULL, prev_nodes_array, _k, g);
        return; // return here since regardless of whether ProcessGraphlet has passed or not, prev_nodes_array is already of size k so we should terminate the recursion
    }

    // Populate next_step with the following algorithm
    SET *next_step; // will contain the set of neighbors of all nodes in prev_nodes_array, excluding nodes actually in prev_nodes_array
    next_step = SetAlloc(G->n);

    // for all nodes in prev_nodes_array...
    for(i=0; i<prev_nodes_count; i++) {
        // loop through all of their neighbors and...
        for(j=0; j<G->degree[prev_nodes_array[i]]; j++) {
            neigh = G->neighbor[prev_nodes_array[i]][j];
            // if the neighbor is not in prev_nodes_array add it to the set
            if(!arrayIn(prev_nodes_array, prev_nodes_count, neigh)) {
                SetAdd(next_step, neigh); // the SET takes care of deduplication
            }
        }
    }

    // create and sort next step heur arr
    int next_step_count = SetCardinality(next_step);
    int next_step_arr[next_step_count];
#if PARANOID_ASSERTS
    assert(SetToArray(next_step_arr, next_step) == next_step_count);
#endif
    SetFree(next_step); // now that we have the next_step_arr, we no longer need the SET next_step

    // populate next_step_nwhn_arr with all nodes in next_step_arr along with their heuristic values provided in heur_arr and their names
    node_whn next_step_nwhn_arr[next_step_count]; // we only need this so that we're able to sort the array
    for (i = 0; i < next_step_count; ++i) {
        int curr_node = next_step_arr[i];
        next_step_nwhn_arr[i].node = curr_node;
        next_step_nwhn_arr[i].heur = heur_arr[curr_node];
        next_step_nwhn_arr[i].name = _nodeNames[curr_node];
    }

    int (*comp_func)(const void*, const void*);

    if (_alphabeticTieBreaking) {
        comp_func = nwhn_des_alph_comp_func;
    } else {
        comp_func = nwhn_des_rev_comp_func;
    }

    qsort((void*)next_step_nwhn_arr, next_step_count, sizeof(node_whn), comp_func); // sort by heuristic first and name second

    // Loop through neighbor nodes with Top N (-lDEGN) distinct heur values
    // If there are multiple nodes with the same heur value (which might happen with degree), we need to expand to all of them because randomly picking one to expand to would break determinism
    // If -lDEGN flag is not given, then will loop through EVERY neighbor nodes in descending order of their degree.
    int num_total_distinct_values = 0;
    double old_heur = -1; // TODO, fix this so that it's not contingent upon heuristics always being >= 0
    i = 0;
    while (i < next_step_count) {
        node_whn next_step_nwhn = next_step_nwhn_arr[i];
        double curr_heur = next_step_nwhn.heur;
        if (curr_heur != old_heur) {
            ++num_total_distinct_values;
        }
        ++i;
    }
    int num_distinct_values_to_skip = (int)(num_total_distinct_values * _topThousandth) / 1000; // algo=base
    // int num_distinct_values_to_skip = _k - prev_nodes_count - 1; // algo=stairs

    int num_distinct_values = 0;
    old_heur = -1; // TODO, fix this so that it's not contingent upon heuristics not being -1
    i = 0;
    while (i < next_step_count) {
        node_whn next_step_nwhn = next_step_nwhn_arr[i];
        double curr_heur = next_step_nwhn.heur;
        if (curr_heur != old_heur) {
            ++num_distinct_values;
        }
        old_heur = curr_heur;

        // continue for the heurs we skip
        if (num_distinct_values <= num_distinct_values_to_skip) {
            ++i;
            continue;
        }

        // break once we've gotten enough distinct heur values
        if (_numWindowRepLimit != 0 && num_distinct_values - num_distinct_values_to_skip > _numWindowRepLimit) {
            break;
        }

        // perform the standard DFS step of set next, recurse with size + 1, and then unset next
        prev_nodes_array[prev_nodes_count] = next_step_nwhn.node;
        SampleGraphletIndexAndPrint(G, prev_nodes_array, prev_nodes_count + 1, heur_arr);
        // the "unset" step is commented out for efficiency since it's not actually necessary, but the comment improves readability
        // prev_nodes_array[prev_nodes_count] = 0;
        ++i;
    }
}

int RunBlantFromGraph(int k, int numSamples, GRAPH *G)
{
    int i,j, windowRepInt, D;
    char perm[MAX_K+1];
    assert(k <= G->n);
    SET *V = SetAlloc(G->n);
    SET *prev_node_set = SetAlloc(G->n);
    SET *intersect_node = SetAlloc(G->n);
    TINY_GRAPH *empty_g = TinyGraphAlloc(k); // allocate it here once, so functions below here don't need to do it repeatedly
    int varraySize = _windowSize > 0 ? _windowSize : MAX_K + 1;
    unsigned Varray[varraySize];
    InitializeConnectedComponents(G);

    if ((_outputMode != indexGraphlets && _outputMode != indexGraphletsRNO && _outputMode != indexOrbits))
	    Fatal("currently only -mi and -mj output modes are supported for INDEX and EDGE_COVER sampling methods");

    int count = 0;
    int prev_nodes_array[_k];

    // Get heuristic values based on orbit number, if ODV file provided
    double heuristicValues[G->n];
    getDoubleDegreeArr(heuristicValues, G); // since heuristic values are doubles, we need to convert degree values to doubles

    int percentToPrint = 1;
    node_whn nwhn_arr[G->n]; // nodes sorted first by the heuristic function and then either alphabetically or reverse alphabetically

    // fill node order array with base values
    for (i = 0; i < G->n; i++) {
        nwhn_arr[i].node = i;
        nwhn_arr[i].heur = heuristicValues[i];
        nwhn_arr[i].name = _nodeNames[i];
    }

    // sort array
    int (*comp_func)(const void*, const void*);

    if (_alphabeticTieBreaking) {
        comp_func = nwhn_des_alph_comp_func;
    } else {
        comp_func = nwhn_des_rev_comp_func;
    }

    qsort((void*)nwhn_arr, G->n, sizeof(node_whn), comp_func);

    for(i=0; i<G->n; i++) {
        prev_nodes_array[0] = nwhn_arr[i].node;

        SampleGraphletIndexAndPrint(G, prev_nodes_array, 1, heuristicValues);
        count = 0;

        if (i * 100 / G->n >= percentToPrint) {
            fprintf(stderr, "%d%% done\n", percentToPrint);
            ++percentToPrint;
        }
    }

#if PARANOID_ASSERTS // no point in freeing this stuff since we're about to exit; it can take significant time for large graphs.
	Free(_graphletDegreeVector[i]);
    TinyGraphFree(empty_g);
#endif
    SetFree(V);
    SetFree(prev_node_set);
    SetFree(intersect_node);
    return 0;
}


static FILE *fpThreads[MAX_POSSIBLE_THREADS]; // these will be the pipes reading output of the parallel blants

void BlantAddEdge(int v1, int v2)
{
    if(!_pairs) _pairs = Malloc(2*_maxEdges*sizeof(_pairs[0]));
    assert(_numEdges <= _maxEdges);
    if(_numEdges >= _maxEdges)
    {
	_maxEdges *=2;
	_pairs = Realloc(_pairs, 2*_maxEdges*sizeof(int));
    }
    _numNodes = MAX(_numNodes, v1+1); // add one since, for example, if we see a node numbered 100, numNodes is 101.
    _numNodes = MAX(_numNodes, v2+1);
    _pairs[2*_numEdges] = v1;
    _pairs[2*_numEdges+1] = v2;
    if(_pairs[2*_numEdges] == _pairs[2*_numEdges+1])
	Fatal("BlantAddEdge: edge %d (%d,%d) has equal nodes; cannot have self-loops\n", _numEdges, v1, v2);
    if(_pairs[2*_numEdges] > _pairs[2*_numEdges+1])
    {
	int tmp = _pairs[2*_numEdges];
	_pairs[2*_numEdges] = _pairs[2*_numEdges+1];
	_pairs[2*_numEdges+1] = tmp;
    }
    assert(_pairs[2*_numEdges] < _pairs[2*_numEdges+1]);
    _numEdges++;
}

const char * const USAGE_SHORT =
"BLANT (Basic Local Alignment of Network Topology): sample graphlets of up to 8 nodes from a graph.\n"\
"USAGE: blant [OPTIONS] -k numNodes -n numSamples graphInputFile\n"\
" Common options: (use -h for longer help)\n"\
"    -s samplingMethod (default MCMC; NBE, EBE, RES, AR, INDEX, EDGE_COVER)\n"\
"    -m{outputMode} (default o=ODV; g=GDV, f=frequency, i=index, r=root(used only for INDEX), d=distribution of neighbors\n"\
"    -d{displayModeForCanonicalIDs} (o=ORCA, j=Jesse, b=binaryAdjMatrix, d=decimal, i=integerOrdinal)\n"\
"    -r seed (integer)\n\n"\
"    -t N[:M]: (CURRENTLY BROKEN): use threading (parallelism); break the task up into N jobs (default 1) allowing\n"\
"        at most M to run at one time; M can be anything from 1 to a compile-time-specified maximum possible value\n"\
"        (MAX_POSSIBLE_THREADS in blant.h), but defaults to 4 to be conservative.";

const char * const USAGE_LONG =
"BLANT: Basic Local Alignment of Network Topology (work in progress)\n"\
"PURPOSE: sample graphlets of up to 8 nodes from a graph. Default output is similar to ORCA, though via stochastic sampling\n"\
"    rather than exaustive enumeration. Our APPROXIMATE results come MUCH faster than ORCA on large or dense networks.\n"\
"USAGE: blant [OPTIONS] -k numNodes -n numSamples graphInputFile\n"\
"where the following are REQUIRED:\n"\
"    numNodes is an integer 3 through 8 inclusive, specifying the size (in nodes) of graphlets to sample;\n"\
"    numSamples is the number of graphlet samples to take (large samples are recommended), except in INDEX sampling mode,\n"\
"	where it specifies the maximum number of samples to take from each node in the graph.\n"\
"    samplingMethod is:\n"\
"	MCMC (default): Markov Chain Monte Carlo: Build the first set S of k nodes using NBE; then randomly remove and add\n"\
"         nodes to S using an MCMC graph walking algorithm with restarts; gives asymptotically correct relative frequencies\n"\
"         when using purely counting modes like -m{o|g|f}, but biased counts in indexing modes like -m{i|j} since we remove\n"\
"           duplicates in indexing modes.)\n"\
"	NBE (Node-Based Expansion): pick a node at random and add it to the node set S; add new nodes by choosing uniformly\n"\
"         at random from all nodes one step outside S. (Fast on sparse networks, slightly biased counts)\n"\
"	EBE (Edge-Based Expansion): pick an edge at random and add its two nodes to S; add nodes to S by picking an edge\n"\
"         uniformly at random from those emanating from S. (faster than NBE on dense networks, but more biased)\n"\
"	RES (Lu Bressan's REServoir sampling): also asymptotically correct but much slower than MCMC.\n"\
"	AR (Accept-Reject): EXTREMELY SLOW but asymptotically correct: pick k nodes entirely at random, reject if\n"\
"	  resulting graphlet is disconnected (vast majority of such grpahlets are disconnected, thus VERY SLOW)\n"\
"	INDEX: deterministic: for each node v in the graph, build a topologically deterministic set of k-graphlets to\n"\
"         be used as indices for seed-and-extend local alignments (using, eg., our onw Dijkstra-inspired local aligner--\n"\
"         see Dijkstra diretory). When using INDEX sampling, the -n command-line option specifies the maximum number\n"\
"         of index entries per starting node v.\n"\
"	EDGE_COVER: starting with E=all edges in the input graph, pick one edge in E and build a k-graphlet using EBE;\n"\
"         then subtract ALL its edges from E. Continue until E is empty. The goal is to output a short list of graphlets\n"\
"         that, in comglomerate, cover each edge in E at least once.\n"\
"    graphInputFile: graph must be in one of the following formats with its extension name:\n"\
"	Edgelist (.el), LEDA(.leda), GML (.gml), GraphML (.xml), LGF(.lgf), CSV(.csv)\n"\
"	(extensions .gz and .xz are automatically decompressed using gunzip and unxz, respectively)\n"\
"	Duplicate edges (either direction) and self-loops should be removed!\n"\
"COMMON OPTIONS:\n"\
"    -m{outputMode}, where {outputMode} is a single character, one of:\n"\
"	o = the default, which is ODV (Orbit Degree Vector), identical to ORCA (commonly though incorrectly called a GDV)\n"\
"	g = GDV (Graphlet Degree Vector) Note this is NOT what is commonly called a GDV, which is actually an ODV (above).\n"\
"	NOTE: the difference is that an ODV counts the number of nodes that touch all possible *orbits*, while a GDV lists\n"\
"		only the smaller vector of how many nodes touch each possible *graphlet* (independent of orbit).\n"\
"	f = graphlet {f}requency, similar to Relative Graphlet Frequency, produces a raw count across our random samples.\n"\
"	    sub-option -mf{freqDispMode} can be i(integer or count) or d(decimal or concentration)\n"\
"	i = {i}ndex: each line is a graphlet with columns: canonical ID, then k nodes in canonical order; useful since\n"\
"	    two lines with the same first column constitutes a PERFECT k-node local alignment between the two graphlets.\n"\
"	r = index with {r}oot node orbit: each line is a canonical ID + the orbit of the root node, then k nodes in canonical order; produces better seeds when the index is queried by the alignment algorithm\n"\
"	d = graphlet neighbor {D}istribution\n"\
"    -d{displayMode} [no default--MANDATORY for indexing modes]: single character controls how canonical IDs are displayed:\n"\
"	o = ORCA numbering\n"\
"	j = JESSE numbering\n"\
"	b = explicit binary representation of the half-adjacency matrix of the canonical graphlet\n"\
"	d = decimal (base-10) integer representation of the above binary\n"\
"	i = integer ordinal = sorting the above integers and numbering them 0, 1, 2, 3, etc.\n"\
"Less Common OPTIONS:\n"\
"    -t N[:M]: use threading (parallelism); break the task up into N jobs (default 1) allowing at most M to run at one time.\n"\
"       M can be anything from 1 to a compile-time-specified maximum possible value (MAX_POSSIBLE_THREADS in blant.h),\n"\
"       but defaults to 4 to be conservative.\n"\
"    -r seed: pick your own random seed\n"\
"    -w windowSize: DEPRECATED. (use '-h' option for more)\n"\
"	-p windowRepSamplingMethod: (DEPRECATED) one of the below\n"\
"	    MIN (Minimizer); MAX (Maximizer); DMIN (Minimizer With Distance); DMAX (Maximizer with Distance);\n"\
"	    LFMIN (Least Frequent Minimizer); LFMAX (Least Frequent Maximizer)\n"\
"	-P windowRepIterationMethods is one of: COMB (Combination) or DFS\n" \
"	-l windowRepLimitMethod is one of: [suffix N: limit to Top N satisfied graphlets]\n"\
"	    DEG (graphlet Total Degree); EDGE (1-step away numEdges)\n"\
"   -M = multiplicity, meaning max allowed number of ambiguous permutations in found graphlets (M=0 is a special case and means no max)\n" \
"   -T = top percent to expand to in -sINDEX sampling method (default 0)\n" \
"   -o = the orbit to use for the heuristic function\n" \
"   -f = the .orca4 file for the network\n" \
"   -a = sets whether ties are broken alphabetically or reverse alphabetically"\
;

// The main program, which handles multiple threads if requested.  We simply fire off a bunch of parallel
// blant *processes* (not threads, but full processes), and simply merge all their outputs together here
// in the parent.
int main(int argc, char *argv[])
{
    int i, j, opt, numSamples=0, multiplicity=1;
    double windowRep_edge_density = 0.0;
    int exitStatus = 0;

    if(argc == 1)
    {
	puts(USAGE_SHORT);
	exit(1);
    }

    _JOBS = 1;
    _MAX_THREADS = 4;

    _k = 0; _k_small = 0;

    int odv_fname_len = 0;

    while((opt = getopt(argc, argv, "hm:d:t:r:s:c:k:K:o:f:e:g:w:p:P:l:n:M:T:a:")) != -1)
    {
	switch(opt)
	{
	case 'h':
	    printf("%s\n", USAGE_LONG);
	    printf("Note: current TSET size is %ld bits\n", 8*sizeof(TSET));
	    exit(1); break;
	case 'm':
	    if(_outputMode != undef) Fatal("tried to define output mode twice");
	    switch(*optarg)
	    {
	    case 'm': _outputMode = indexMotifs; break;
	    case 'M': _outputMode = indexMotifOrbits; break;
	    case 'i': _outputMode = indexGraphlets; break;
	    case 'r': _outputMode = indexGraphletsRNO; break;
	    case 'j': _outputMode = indexOrbits; break;
	    case 'f': _outputMode = graphletFrequency;
		switch (*(optarg + 1))
		{
		    case 'i': _freqDisplayMode = count; break;
		    case 'd': _freqDisplayMode = concentration; break;
		    case '\0': _freqDisplayMode = freq_display_mode_undef; break;
		    default: Fatal("-mf%c: unknown frequency display mode;\n"
		    "\tmodes are i=integer(count), d=decimal(concentration)", *(optarg + 1));
		    break;
		}
	    break;
	    case 'g': _outputMode = outputGDV; break;
	    case 'o': _outputMode = outputODV; break;
	    case 'd': _outputMode = graphletDistribution; break;
	    default: Fatal("-m%c: unknown output mode \"%c\"", *optarg,*optarg);
	    break;
	    }
	    break;
	case 'd':
	    if (_displayMode != undefined) Fatal("tried to define canonical display mode twice");
	    switch(*optarg)
	    {
	    case 'b': _displayMode = binary; break;
	    case 'd': _displayMode = decimal; break;
	    case 'i': _displayMode = ordinal; break;
	    case 'j': _displayMode = jesse; break;
	    case 'o': _displayMode = orca; break;
	    default: Fatal("-d%c: unknown canonical display mode:n"
		    "\tmodes are i=integer ordinal, d=decimal, b=binary, o=orca, j=jesse", *optarg);
	    break;
	    }
	    break;
	case 't': assert(1==sscanf(optarg, "%d", &_JOBS));
	    _MAX_THREADS = _JOBS;
	    assert(1 <= _JOBS && _MAX_THREADS <= MAX_POSSIBLE_THREADS);
	    break;
	case 'r': _seed = atoi(optarg); if(_seed==-1)Apology("seed -1 ('-r -1' is reserved to mean 'uninitialized'");
	    break;
	case 'k': _k = atoi(optarg);
	    if (!(3 <= _k && _k <= 8)) Fatal("%s\nERROR: k [%d] must be between 3 and 8\n%s", USAGE_SHORT, _k);
	    break;
	case 'w': _window = true; _windowSize = atoi(optarg); break;
	case 'p':
	    if (_windowSampleMethod != -1) Fatal("Tried to define window sampling method twice");
	    else if (strncmp(optarg, "DMIN", 4) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_MIN_D;
	    else if (strncmp(optarg, "DMAX", 4) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_MAX_D;
	    else if (strncmp(optarg, "MIN", 3) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_MIN;
	    else if (strncmp(optarg, "MAX", 3) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_MAX;
	    else if (strncmp(optarg, "LFMIN", 5) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_LEAST_FREQ_MIN;
	    else if (strncmp(optarg, "LFMAX", 5) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_LEAST_FREQ_MAX;
	    else if (strncmp(optarg, "DEGMAX", 6) == 0)
		_windowSampleMethod = WINDOW_SAMPLE_DEG_MAX;
	    else
		Fatal("Unrecognized window searching method specified. Options are: -p[u|U]{MIN|MAX|DMIN|DMAX|LFMIN|LFMAX|DEGMAX}\n");
	    break;
	case 'P':
		if (strncmp(optarg, "COMB", 4) == 0)
			_windowIterationMethod = WINDOW_ITER_COMB;
		else if (strncmp(optarg, "DFS", 3) == 0)
			_windowIterationMethod = WINDOW_ITER_DFS;
		else
			Fatal("Unrecognized window Iteration method specified. Options are: -P{COMB|DFS}\n");
		break;
	case 'l':
		if (_windowRep_limit_method != WINDOW_LIMIT_UNDEF) Fatal("Tried to define window limiting method twice");
		if (strncmp(optarg, "n", 1) == 0 || strncmp(optarg, "N", 1) == 0) {
		    _windowRep_limit_neglect_trivial = true; optarg += 1;
		}
		if (strncmp(optarg, "DEG", 3) == 0) {
			_windowRep_limit_method = WINDOW_LIMIT_DEGREE; optarg += 3;
		}
		else if (strncmp(optarg, "EDGE", 4) == 0) {
			_windowRep_limit_method = WINDOW_LIMIT_EDGES; optarg += 4;
		}
		else
			Fatal("Unrecognized window limiting method specified. Options are: -l{DEG}{EDGE}{limit_num}\n");
		_numWindowRepLimit = atoi(optarg);
		if (!_numWindowRepLimit) {_numWindowRepLimit = 10; _numWindowRepArrSize = _numWindowRepLimit;}
		_windowRep_limit_heap = HeapAlloc(_numWindowRepLimit, asccompFunc, NULL);
		break;
	case 'n': numSamples = atoi(optarg);
		char lastChar = optarg[strlen(optarg)-1];
		if(!isdigit(lastChar))
		    switch(lastChar) {
		    case 'b': case 'B': case 'g': case 'G': numSamples *= 1024; // do NOT break, fall through
		    case 'm': case 'M': numSamples *= 1024;
		    case 'k': case 'K': numSamples *= 1024; break;
		    default: Fatal("%s\nERROR: numSamples can be appended by k, m, b, or g but not %c\n%s", USAGE_SHORT, lastChar);
		    break;
		    }
		if(numSamples < 0) Fatal("%s\nFatal Error: numSamples [%d] must be non-negative", USAGE_SHORT, numSamples);
		//fprintf(stderr, "numSamples set to %d\n", numSamples);
		break;
	case 'M': multiplicity = atoi(optarg);
	    if(multiplicity < 0) Fatal("%s\nERROR: multiplicity [%d] must be non-negative\n", USAGE_SHORT, multiplicity);
    case 'T': _topThousandth = atoi(optarg);
	    break;
	case 'o':
        _orbitNumber = atoi(optarg);
        break;
    case 'f':
        odv_fname_len = strlen(optarg);
        _odvFile = malloc(sizeof(char) * odv_fname_len);
        strncpy(_odvFile, optarg, odv_fname_len);
        break;
    case 'a':
        _alphabeticTieBreaking = atoi(optarg) != 0;
        break;
	default: Fatal("%s\nERROR: unknown option %c", USAGE_SHORT, opt);
    }
    }

    if (_k <= 5) Fatal("k is %d but must be at least 6 because there are no unambiguous graphlets for k<=5",_k); // PATNOTE: keep me

    if(_seed == -1) _seed = GetFancySeed(false);
    // This only seeds the main thread; sub-threads, if they exist, are seeded later by "stealing"
    // exactly _THREADS-1 values from near the beginning of this main random stream.
    RandomSeed(_seed);

    if(_outputMode == undef) _outputMode = outputODV; // default to the same thing ORCA and Jesse us
	if (_freqDisplayMode == freq_display_mode_undef) // Default to integer(count)
		_freqDisplayMode = count;

    FILE *fpGraph;
    int piped = 0;
    if(!argv[optind])
    {
	fpGraph = stdin;
	if(isatty(0)) Warning("reading graph input file from terminal, press ^D to finish");
    }
    else {
	char *graphFileName = argv[optind];
	fpGraph = readFile(graphFileName, &piped);
	if(!fpGraph) Fatal("cannot open graph input file '%s'\n", argv[optind]);
	optind++;
    }

    SetBlantDir(); // Needs to be done before reading any files in BLANT directory
    SetGlobalCanonMaps(); // needs _k to be set
    LoadMagicTable(); // needs _k to be set

    // This section of code computes the info necessary to implement "multiplicity" mode ('M' above), which controls
    // whether to output a sampled graphlet at all during INDEXING based on how much "ambiguity" there is in its
    // local alignment. The more permutations of the graphlet there are, the less useful it is for seeding local
    // alignments and therefore the less useful as a database index entry.
    _windowRep_allowed_ambig_set = SetAlloc(_numCanon);
    SET *orbit_temp = SetAlloc(_numOrbits);
    for(i=0; i<_numCanon; i++) {
        if (SetIn(_connectedCanonicals, i)) {
            // calculate number of permutations for the given canonical graphlet
	    // (loop through all unique orbits and count how many times they appear)
            // the formula is for every unique orbit, multiply the number of permutations by
	    // the factorial of how many appearances that unique orbit has
            // if there is one orbit with three nodes and a second orbit with 2 nodes,
	    // the number of permutations would be (3!)(2!)
            // NOTE: this is not the most efficient algorithm since it doesn't use hash tables.
	    // I didn't want to overcomplicate it because it only happens once per run.
            // however, if speed is important (this currently takes about 5 seconds on k=8) this can be sped up
            for(j=0; j<_k; j++) SetAdd(orbit_temp, _orbitList[i][j]);
            unsigned uniq_orbits[_k];
            unsigned num_uniq_orbits = SetToArray(uniq_orbits, orbit_temp);
            unsigned uniq_orbit_i;
            unsigned total_orbit_perms = 1;

            for (uniq_orbit_i=0; uniq_orbit_i<num_uniq_orbits; uniq_orbit_i++) {
                unsigned orbit_appearances = 0;
                unsigned orbit_i = 0;

                for (orbit_i=0; orbit_i<_k; orbit_i++) {
                    if (_orbitList[i][orbit_i] == uniq_orbits[uniq_orbit_i]) {
                        orbit_appearances++;
                        total_orbit_perms *= orbit_appearances;
                    }
                }
            }

            // I know it's inefficient to put multiplicity here instead of around the whole orbit perm calculation code but
	    // it increases readability, at least until orbit perm calculation is put into a function
            if(multiplicity == 0 || total_orbit_perms <= multiplicity) { // multiplicity = 0 means any ambiguity is allowed
                SetAdd(_windowRep_allowed_ambig_set, i);
            }
            SetEmpty(orbit_temp);
        }
    }
    SetFree(orbit_temp);

    // Read network using native Graph routine.
    GRAPH *G = GraphReadEdgeList(fpGraph, SPARSE, _supportNodeNames);
    if(_supportNodeNames)
    {
	assert(G->name);
	_nodeNames = G->name;
    }
    if(fpGraph != stdin) closeFile(fpGraph, &piped);

	exitStatus = RunBlantFromGraph(_k, numSamples, G);
    GraphFree(G);
    return exitStatus;
}
