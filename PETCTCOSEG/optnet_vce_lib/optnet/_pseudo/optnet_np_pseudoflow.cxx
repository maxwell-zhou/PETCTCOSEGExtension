/*
 ==========================================================================
 |   
 |   $Id: optnet_np_pseudoflow.cxx 21 2008-01-12 15:52:31Z Qi Song $
 |
 |   Written by Qi Song <qi-song@uiowa.edu>
 |   Department of Electrical and Computer Engineering
 |   University of Iowa
 |   
 ==========================================================================
 */

#ifndef ___OPTNET_PSEUDOFLOW_CXX___
#   define ___OPTNET_PSEUDOFLOW_CXX___

#   include <optnet/_pseudo/optnet_np_pseudoflow.hxx>
#   include <limits>

#   ifdef max       // The max macro may interfere with
#       undef max   //   std::numeric_limits::max().
#   endif  
#   define MAX_VALUE 10000000000
//

namespace optnet{

///////////////////////////////////////////////////////////////////////////
template <typename _Cap>
optnet_pseudoflow<_Cap>::optnet_pseudoflow()
{
    numNodes=0;
    numArc1s=0;
	arcIndex=0;
    source = 1;
    sink = 2;
    numParams = 0;

    highestStrongLabel = 1;

	labelCount = NULL;
    adjacencyList = NULL;
    strongRoots = NULL;

	labelList = NULL;  
     numPushes = 0;
     numMergers = 0;
     numRelabels = 0;
     numGaps = 0;
     numArc1Scans = 0;
	 m_lambda=1;

}


///////////////////////////////////////////////////////////////////////////
template <typename _Cap>
bool
optnet_pseudoflow<_Cap>::create(size_type size_x, size_type size_y, size_type size_z, size_type size_s)
                                   
{   size_type i;
	m_x = size_x;
	m_y = size_y;
	m_z = size_z;
	m_s = size_s;
	m_numcols = size_x * size_y * size_s;
	numNodes = m_z * m_numcols + 2;
	numParams = 1;    
    adjacencyList = new Node [numNodes];
	strongRoots = new Root [numNodes];
	labelCount = new size_type [numNodes];
	labelList = new size_type [numNodes];

	for ( i = 0; i < numNodes; ++i)
	{
		initializeRoot (&strongRoots[i]);
		initializeNode (&adjacencyList[i], (i+1));
		labelCount[i] = 0;
		labelList[i]=0;
	}

	return true;

}
/////////////////////////////////
template <typename _Cap>
bool
optnet_pseudoflow<_Cap>::create(size_type numpc, size_type numcols)
                                   
{   size_type i;
	m_colsize = numpc;
	m_numcols = numcols;
	numNodes = m_colsize * m_numcols + 2;
	numParams = 1;    
    adjacencyList = new Node [numNodes];
	strongRoots = new Root [numNodes];
	labelCount = new size_type [numNodes];
	labelList = new size_type [numNodes];

	for ( i = 0; i < numNodes; ++i)
	{
		initializeRoot (&strongRoots[i]);
		initializeNode (&adjacencyList[i], (i+1));
		labelCount[i] = 0;
		labelList[i] = 0;
	}

	return true;

}
//////////////////////////


template <typename _Cap>
void optnet_pseudoflow<_Cap>::add_arc(size_type tail_x, size_type tail_y, size_type tail_z, size_type head_x, size_type head_y, size_type head_z)
{
	//int i,j;
	size_type from, to;
    Arc1 *ac=new Arc1;
	initializeArc1 (ac);
	from= (tail_x*m_y+tail_y)*m_z+tail_z+3;
	to=(head_x*m_y+head_y)*m_z+head_z+3;

	ac->from = &adjacencyList[from-1];
	ac->to = &adjacencyList[to-1];

	if ((from!= source) && (to!=sink))
	{
		ac->capacities= new capacity_type;
		ac->capacities[0]=MAX_VALUE;
		ac->capacity= ac->capacities[0];
	}
    ++ arcIndex;
	++ ac->from->numAdjacent;
	++ ac->to->numAdjacent;
	Arc1List.push_back(ac);
}
template <typename _Cap>
void optnet_pseudoflow<_Cap>::add_arc(size_type tail_x, size_type tail_y, size_type tail_z, size_type tail_s, size_type head_x, size_type head_y, size_type head_z, size_type head_s)
{
	//int i,j;
	size_type from, to;
    Arc1 *ac=new Arc1;
	initializeArc1 (ac);
	from= m_x*m_y*m_z*tail_s+(tail_x*m_y+tail_y)*m_z+tail_z+3;
	to= m_x*m_y*m_z*head_s+(head_x*m_y+head_y)*m_z+head_z+3;

	ac->from = &adjacencyList[from-1];
	ac->to = &adjacencyList[to-1];

	if ((from!= source) && (to!=sink))
	{
		ac->capacities= new capacity_type;
		ac->capacities[0]=MAX_VALUE;
		ac->capacity= ac->capacities[0];
	}
    ++ arcIndex;
	++ ac->from->numAdjacent;
	++ ac->to->numAdjacent;
	Arc1List.push_back(ac);
}

////////////////////////////////
template <typename _Cap>
void optnet_pseudoflow<_Cap>::add_arc(size_type tail_index_npc, size_type tail_index_col, size_type head_index_npc, size_type head_index_col)
{
	//int i,j;
	size_type from, to;
    Arc1 *ac=new Arc1;
	initializeArc1 (ac);
	from= tail_index_col*m_colsize+tail_index_npc+3;
	to=head_index_col*m_colsize+head_index_npc+3;

	ac->from = &adjacencyList[from-1];
	ac->to = &adjacencyList[to-1];

	if ((from!= source) && (to!=sink))
	{
		ac->capacities= new capacity_type;
		ac->capacities[0]=MAX_VALUE;
		ac->capacity= ac->capacities[0];
	}
    ++ arcIndex;
	++ ac->from->numAdjacent;
	++ ac->to->numAdjacent;
	Arc1List.push_back(ac);
}

/////////////////////////////////
template <typename _Cap>
void optnet_pseudoflow<_Cap>::add_st_arc(capacity_type s, capacity_type t, size_type index_x, size_type index_y, size_type index_z, size_type index_s)

{
	size_type from, to;
    

	if ( s !=0 )
	{ from = source;
	  to = m_x * m_y * m_z * index_s + ( index_x * m_y + index_y ) * m_z + index_z + 3;
	  Arc1 *ac = new Arc1;
	  initializeArc1 (ac);
	  ac->from = &adjacencyList[ from - 1 ];
	  ac->to = &adjacencyList[ to - 1 ];
	  ac->capacities = new capacity_type [ numParams ];
	  ac->capacities[ 0 ] = s;
	  ac->capacity = ac->capacities[ 0 ];
	  ++ arcIndex;
      ++ ac->from->numAdjacent;
	  ++ ac->to->numAdjacent;
	  Arc1List.push_back( ac );	
	}
	
	if ( t != 0 )
	{
	  from = m_x * m_y * m_z * index_s + ( index_x * m_y + index_y ) * m_z + index_z + 3;
	  to = sink;
	  Arc1 *ac = new Arc1;
	  initializeArc1 (ac);
	  ac->from = &adjacencyList[ from - 1 ];
	  ac->to = &adjacencyList[ to - 1 ];
	  ac->capacities = new capacity_type [ numParams ];
	  ac->capacities[ 0 ] = t;
	  ac->capacity = ac->capacities[ 0 ];
	  ++ arcIndex;
      ++ ac->from->numAdjacent;
	  ++ ac->to->numAdjacent;
	  Arc1List.push_back( ac );	
	}

		
}
//////////////////////////////////
template <typename _Cap>
void optnet_pseudoflow<_Cap>::add_st_arc(capacity_type s, capacity_type t, size_type index_npc, size_type index_col)

{
	size_type from, to;


	if ( s != 0 )
	{ from = source;
	  to = index_col * m_colsize + index_npc + 3;
	  Arc1 *ac=new Arc1;
	  initializeArc1 (ac);
	  ac->from = &adjacencyList[ from - 1 ];
	  ac->to = &adjacencyList[ to - 1 ];
	  ac->capacities= new capacity_type [ numParams ];
	  ac->capacities[ 0 ] = s;
	  ac->capacity = ac->capacities[ 0 ];
	  ++ arcIndex;
      ++ ac->from->numAdjacent;
	  ++ ac->to->numAdjacent;
	  Arc1List.push_back( ac );	
	}

	if ( t != 0 )
	{
	  from = index_col * m_colsize + index_npc + 3;
	  to = sink;
	  Arc1 *ac=new Arc1;
	  initializeArc1 (ac);
	  ac->from = &adjacencyList[ from - 1 ];
	  ac->to = &adjacencyList[ to - 1 ];
	  ac->capacities = new capacity_type [ numParams ];
	  ac->capacities[ 0 ] = t;
	  ac->capacity = ac->capacities[ 0 ];
	  ++ arcIndex;
      ++ ac->from->numAdjacent;
	  ++ ac->to->numAdjacent;
	  Arc1List.push_back( ac );	

	}

		
}
///////////////////////////////////
template <typename _Cap>
void optnet_pseudoflow<_Cap>:: add_arc_cost(capacity_type edge_cost,size_type tail_index_npc, size_type tail_index_col, size_type head_index_npc, size_type head_index_col)
{
	size_type from, to;
    Arc1 *ac=new Arc1;
	initializeArc1 (ac);
	from= tail_index_col*m_colsize+tail_index_npc+3;
	to=head_index_col*m_colsize+head_index_npc+3;

	ac->from = &adjacencyList[from-1];
	ac->to = &adjacencyList[to-1];

	if ((from!= source) && (to!=sink))
	{
		ac->capacities= new capacity_type;
		ac->capacities[0]=edge_cost;
		ac->capacity= ac->capacities[0];
	}
    ++ arcIndex;
	++ ac->from->numAdjacent;
	++ ac->to->numAdjacent;
	Arc1List.push_back(ac);
}
///////////////////////////////////
template <typename _Cap>
void optnet_pseudoflow<_Cap>::add_arc_cost(capacity_type edge_cost,size_type tail_x, size_type tail_y, size_type tail_z, size_type head_x, size_type head_y, size_type head_z)
{
	size_type from, to;
    Arc1 *ac=new Arc1;
	initializeArc1 (ac);
	from= (tail_x*m_y+tail_y)*m_z+tail_z+3;
	to=(head_x*m_y+head_y)*m_z+head_z+3;;

	ac->from = &adjacencyList[from-1];
	ac->to = &adjacencyList[to-1];

	if ((from!= source) && (to!=sink))
	{
		ac->capacities= new capacity_type;
		ac->capacities[0]=edge_cost;
		ac->capacity= ac->capacities[0];
	}
    ++ arcIndex;
	++ ac->from->numAdjacent;
	++ ac->to->numAdjacent;
	Arc1List.push_back(ac);
}
///////////////////////////////////
template <typename _Cap>
void optnet_pseudoflow<_Cap>::add_arc_cost(capacity_type edge_cost,size_type tail_x, size_type tail_y, size_type tail_z, size_type tail_s, size_type head_x, size_type head_y, size_type head_z, size_type head_s)
{
	size_type from, to;
    Arc1 *ac=new Arc1;
	initializeArc1 (ac);
	from= (m_x*m_y*m_z)*tail_s+(tail_x*m_y+tail_y)*m_z+tail_z+3;
	to=(m_x*m_y*m_z)*head_s+(head_x*m_y+head_y)*m_z+head_z+3;;

	ac->from = &adjacencyList[from-1];
	ac->to = &adjacencyList[to-1];

	if ((from!= source) && (to!=sink))
	{
		ac->capacities= new capacity_type;
		ac->capacities[0]=edge_cost;
		ac->capacity= ac->capacities[0];
	}
    ++ arcIndex;
	++ ac->from->numAdjacent;
	++ ac->to->numAdjacent;
	Arc1List.push_back(ac);
}
/////////////////////////////////
template <typename _Cap>
void optnet_pseudoflow<_Cap>::prepareList()
{
	size_type i,from, to,capacity;
	for (i=0; i<numNodes; ++i) 
	{
		createOutOfTree (&adjacencyList[i]);
	}
	numArc1s=Arc1List.size();

	for (i=0; i<numArc1s; i++) 
	{
		to = Arc1List[i]->to->number;
		from = Arc1List[i]->from->number;
		capacity = Arc1List[i]->capacity;
		//printf("to is %d\n", to);
		//printf("from is %d\n", from);
        //printf ("capacity is     : %d\n", capacity);
		if (!((source == to) || (sink == from) || (from == to))) 
		{
			if ((source == from) && (to == sink)) 
			{
				Arc1List[i]->flow = capacity;
			}
			else if (from == source)
			{
				addOutOfTreeNode (&adjacencyList[from-1], Arc1List[i]);
			}
			else if (to == sink)
			{
				addOutOfTreeNode (&adjacencyList[to-1], Arc1List[i]);
			}
			else
			{
				addOutOfTreeNode (&adjacencyList[from-1], Arc1List[i]);
			}
		}
	}

}
template <typename _Cap>
 bool optnet_pseudoflow<_Cap>::in_source_set(size_type index_x, size_type index_y, size_type index_z, size_type index_s)
 {
	 int n_index;
	 n_index= m_x*m_y*m_z*index_s+(index_x*m_y+index_y)*m_z+index_z+3;
	 if (labelList[n_index-1]==1)
		 return true;
	 else
		 return false;
 }
///////////////
 template <typename _Cap>
 bool optnet_pseudoflow<_Cap>::in_source_set(size_type index_npc, size_type index_col)
 {
	 int n_index;
	 n_index= index_col*m_colsize+index_npc+3;
	 if (labelList[n_index-1]==1)
		 return true;
	 else
		 return false;
 }
 ///////////////////////////////////////////////////////////////////////////
template <typename _Cap>
void
optnet_pseudoflow<_Cap>::initialize()
{
    printf ("c Pseudoflow algorithm for parametric min cut (version 1.0)\n");
	//readDimacsFileCreateList ();
	prepareList();                       //set as public
	simpleInitialization ();             //set as public

	printf ("c Finished initialization.\n"); fflush (stdout);
	
}



///////////////////////////////////////////////////////////////////////////
template <typename _Cap>
typename optnet_pseudoflow<_Cap>::capacity_type
optnet_pseudoflow<_Cap>::solve()
{
    printf ("c Pseudoflow algorithm for parametric min cut (version 1.0)\n");
	//readDimacsFileCreateList ();
	prepareList();    
	printf ("c Finished list preparing.\n");//set as public
	simpleInitialization ();             //set as public

	printf ("c Finished initialization.\n");
	pseudoflowPhase1 ();

	printf ("c Finished phase 1.\n"); fflush (stdout);

//-----------------------------------------------
	/*BySq
	recoverFlow();
	checkOptimality ();
	*/

	printf ("c Number of nodes     : %d\n", numNodes);
	printf ("c Number of Arc1s      : %ld\n", numArc1s);

	printf ("c Number of Arc1 scans : %lld\n", numArc1Scans);
	printf ("c Number of mergers   : %d\n", numMergers);
	printf ("c Number of pushes    : %lld\n", numPushes);
	printf ("c Number of relabels  : %d\n", numRelabels);
	printf ("c Number of gaps      : %d\n", numGaps);



	//displayBreakpoints ();
	//freeMemory ();  //set as public
	return m_flow;


}

///////////////////////////////////////////////////////////////////////////
//----------------------solve the problem
template <typename _Cap>
 void optnet_pseudoflow<_Cap>::initializeNode (Node *nd, size_type n)
{
	nd->label = 0;
	nd->excess = 0;
	nd->parent = NULL;
	nd->childList = NULL;
	nd->nextScan = NULL;
	nd->nextArc1 = 0;
	nd->numOutOfTree = 0;
	nd->Arc1ToParent = NULL;
	nd->next = NULL;
	nd->prev = NULL;
	nd->visited = 0;
	nd->numAdjacent = 0;
	nd->number = n;
	nd->outOfTree = NULL;
	nd->breakpoint = (numParams+1);
}
template <typename _Cap>
void optnet_pseudoflow<_Cap>::initializeRoot (Root *rt) 
{
	//rt->start = (Node *) malloc (sizeof(Node));
	rt->start = new Node;
	//rt->end = (Node *) malloc (sizeof(Node));
	rt->end= new Node;

	/*
	if ((rt->start == NULL) || (rt->end == NULL))
	{
		printf ("%s Line %d: Out of memory\n", __FILE__, __LINE__);
		exit (1);
	}
	*/

	initializeNode (rt->start, 0);
	initializeNode (rt->end, 0);

	rt->start->next = rt->end;
	rt->end->prev = rt->start;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::freeRoot (Root *rt) 
{
	free(rt->start);
	rt->start = NULL;

	free(rt->end);
	rt->end = NULL;
}
template <typename _Cap>
void optnet_pseudoflow<_Cap>::liftAll (Node *rootNode, size_type theparam) 
{
	Node *temp, *current=rootNode;

	current->nextScan = current->childList;

	-- labelCount[current->label];
	current->label = numNodes;	
	current->breakpoint = (theparam+1);

	for ( ; (current); current = current->parent)
	{
		while (current->nextScan) 
		{
			temp = current->nextScan;
			current->nextScan = current->nextScan->next;
			current = temp;
			current->nextScan = current->childList;

			-- labelCount[current->label];
			current->label = numNodes;
			current->breakpoint = (theparam+1);	
		}
	}
}
template <typename _Cap>
void optnet_pseudoflow<_Cap>::addToStrongBucket (Node *newRoot, Node *rootEnd) 
{
	newRoot->next = rootEnd;
	newRoot->prev = rootEnd->prev;
	rootEnd->prev = newRoot;
	newRoot->prev->next = newRoot;
}
template <typename _Cap>
void optnet_pseudoflow<_Cap>::createOutOfTree (Node *nd)
{
	if (nd->numAdjacent)
	{
		/*
		if ((nd->outOfTree = (Arc1 **) malloc (nd->numAdjacent * sizeof (Arc1 *))) == NULL)
		{
			printf ("%s Line %d: Out of memory\n", __FILE__, __LINE__);
			//exit (1);
		}
		*/
		nd->outOfTree= new Arc1* [nd->numAdjacent];
		
	}
}
template <typename _Cap>
void optnet_pseudoflow<_Cap>::initializeArc1 (Arc1 *ac)
{
	//int i;

	ac->from = NULL;
	ac->to = NULL;
	ac->capacity = 0;
	ac->flow = 0;
	ac->direction = 1;
	ac->capacities = NULL;
	ac->ispara=false;
}
template <typename _Cap>
void optnet_pseudoflow<_Cap>::addOutOfTreeNode (Node *n, Arc1 *out) 
{
	n->outOfTree[n->numOutOfTree] = out;
	++ n->numOutOfTree;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::simpleInitialization (void) 
{
	size_type i, size;
	Arc1 *tempArc1;

	size = adjacencyList[source-1].numOutOfTree;
	for (i=0; i<size; ++i) 
	{
		tempArc1 = adjacencyList[source-1].outOfTree[i];
		tempArc1->flow = tempArc1->capacity;
		tempArc1->to->excess += tempArc1->capacity;
	}

	size = adjacencyList[sink-1].numOutOfTree;
	for (i=0; i<size; ++i)
	{
		tempArc1 = adjacencyList[sink-1].outOfTree[i];
		tempArc1->flow = tempArc1->capacity;
		tempArc1->from->excess -= tempArc1->capacity;
	}

	adjacencyList[source-1].excess = 0;
	adjacencyList[sink-1].excess = 0;

	for (i=0; i<numNodes; ++i) 
	{
		if (adjacencyList[i].excess > 0) 
		{
		    adjacencyList[i].label = 1;
			++ labelCount[1];

			addToStrongBucket (&adjacencyList[i], strongRoots[1].end);
		}
	}

	adjacencyList[source-1].label = numNodes;
	adjacencyList[source-1].breakpoint = 0;
	adjacencyList[sink-1].label = 0;
	adjacencyList[sink-1].breakpoint = (numParams+2);
	labelCount[0] = (numNodes - 2) - labelCount[1];
}

template <typename _Cap>
int optnet_pseudoflow<_Cap>::addRelationship (Node *newParent, Node *child) 
{
	child->parent = newParent;
	child->next = newParent->childList;
	newParent->childList = child;

	return 0;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::breakRelationship (Node *oldParent, Node *child) 
{
	Node *current;

	child->parent = NULL;

	if (oldParent->childList == child) 
	{
		oldParent->childList = child->next;
		child->next = NULL;
		return;
	}

	for (current = oldParent->childList; (current->next != child); current = current->next);

	current->next = child->next;
	child->next = NULL;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::merge (Node *parent, Node *child, Arc1 *newArc1) 
{
	Arc1 *oldArc1;
	Node *current = child, *oldParent, *newParent = parent;

#ifdef STATS
	++ numMergers;
#endif

	while (current->parent) 
	{
		oldArc1 = current->Arc1ToParent;
		current->Arc1ToParent = newArc1;
		oldParent = current->parent;
		breakRelationship (oldParent, current);
		addRelationship (newParent, current);
		newParent = current;
		current = oldParent;
		newArc1 = oldArc1;
		newArc1->direction = 1 - newArc1->direction;
	}

	current->Arc1ToParent = newArc1;
	addRelationship (newParent, current);
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::pushUpward (Arc1 *currentArc1, Node *child, Node *parent, const int resCap) 
{
#ifdef STATS
	++ numPushes;
#endif

	if (resCap >= child->excess) 
	{
		parent->excess += child->excess;
		currentArc1->flow += child->excess;
		child->excess = 0;
		return;
	}

	currentArc1->direction = 0;
	parent->excess += resCap;
	child->excess -= resCap;
	currentArc1->flow = currentArc1->capacity;
	parent->outOfTree[parent->numOutOfTree] = currentArc1;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);

	addToStrongBucket (child, strongRoots[child->label].end);
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::pushDownward (Arc1 *currentArc1, Node *child, Node *parent, int flow) 
{
#ifdef STATS
	++ numPushes;
#endif

	if (flow >= child->excess) 
	{
		parent->excess += child->excess;
		currentArc1->flow -= child->excess;
		child->excess = 0;
		return;
	}

	currentArc1->direction = 1;
	child->excess -= flow;
	parent->excess += flow;
	currentArc1->flow = 0;
	parent->outOfTree[parent->numOutOfTree] = currentArc1;
	++ parent->numOutOfTree;
	breakRelationship (parent, child);

	addToStrongBucket (child, strongRoots[child->label].end);
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::pushExcess (Node *strongRoot) 
{
	Node *current, *parent;
	Arc1 *Arc1ToParent;

	for (current = strongRoot; (current->excess && current->parent); current = parent) 
	{
		parent = current->parent;
		Arc1ToParent = current->Arc1ToParent;
		if (Arc1ToParent->direction)
		{
			pushUpward (Arc1ToParent, current, parent, (Arc1ToParent->capacity - Arc1ToParent->flow)); 
		}
		else
		{
			pushDownward (Arc1ToParent, current, parent, Arc1ToParent->flow); 
		}
	}

	if (current->excess > 0) 
	{
		if (!current->next)
		{
			addToStrongBucket (current, strongRoots[current->label].end);
		}
	}
}

template <typename _Cap>
typename optnet_pseudoflow<_Cap>::Arc1 *
optnet_pseudoflow<_Cap>::findWeakNode (Node *strongNode, Node **weakNode) 
{
	int i, size;
	Arc1 *out;

	size = strongNode->numOutOfTree;

	for (i=strongNode->nextArc1; i<size; ++i) 
	{


		++ numArc1Scans;


		if (strongNode->outOfTree[i]->to->label == (highestStrongLabel-1)) 
		{
			strongNode->nextArc1 = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->to;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		}
		else if (strongNode->outOfTree[i]->from->label == (highestStrongLabel-1)) 
		{
			strongNode->nextArc1 = i;
			out = strongNode->outOfTree[i];
			(*weakNode) = out->from;
			-- strongNode->numOutOfTree;
			strongNode->outOfTree[i] = strongNode->outOfTree[strongNode->numOutOfTree];
			return (out);
		}
	}

	strongNode->nextArc1 = strongNode->numOutOfTree;

	return NULL;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::checkChildren (Node *curNode) 
{
	for ( ; (curNode->nextScan); curNode->nextScan = curNode->nextScan->next)
	{
		if (curNode->nextScan->label == curNode->label)
		{
			return;
		}
		
	}	

	-- labelCount[curNode->label];
	++	curNode->label;
	++ labelCount[curNode->label];
	++ numRelabels;
	curNode->nextArc1 = 0;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::processRoot (Node *strongRoot) 
{
	Node *temp, *strongNode = strongRoot, *weakNode;
	Arc1 *out;

	strongRoot->nextScan = strongRoot->childList;

	if ((out = findWeakNode (strongRoot, &weakNode)))
	{
		merge (weakNode, strongNode, out);
		pushExcess (strongRoot);
		return;
	}

	checkChildren (strongRoot);
	
	while (strongNode)
	{
		while (strongNode->nextScan) 
		{
			temp = strongNode->nextScan;
			strongNode->nextScan = strongNode->nextScan->next;
			strongNode = temp;
			strongNode->nextScan = strongNode->childList;

			if ((out = findWeakNode (strongNode, &weakNode)))
			{
				merge (weakNode, strongNode, out);
				pushExcess (strongRoot);
				return;
			}

			checkChildren (strongNode);
		}

		if ((strongNode = strongNode->parent))
		{
			checkChildren (strongNode);
		}
	}

	addToStrongBucket (strongRoot, strongRoots[strongRoot->label].end);

	++ highestStrongLabel;
}

template <typename _Cap>
typename optnet_pseudoflow<_Cap>::Node *
optnet_pseudoflow<_Cap>::getHighestStrongRoot (const int theparam) 
{
	int i;
	Node *strongRoot;

	for (i=highestStrongLabel; i>0; --i) 
	{
		if (strongRoots[i].start->next != strongRoots[i].end)  
		{
			highestStrongLabel = i;
			if (labelCount[i-1]) 
			{
				strongRoot = strongRoots[i].start->next;
				strongRoot->next->prev = strongRoot->prev;
				strongRoot->prev->next = strongRoot->next;
				strongRoot->next = NULL;
				return strongRoot;				
			}

			while (strongRoots[i].start->next != strongRoots[i].end) 
			{


				++ numGaps;

				strongRoot = strongRoots[i].start->next;
				strongRoot->next->prev = strongRoot->prev;
				strongRoot->prev->next = strongRoot->next;
				liftAll (strongRoot, theparam);
			}
		}
	}

	if (strongRoots[0].start->next == strongRoots[0].end) 
	{
		return NULL;
	}

	while (strongRoots[0].start->next != strongRoots[0].end) 
	{
		strongRoot = strongRoots[0].start->next;
		strongRoot->next->prev = strongRoot->prev;
		strongRoot->prev->next = strongRoot->next;

		strongRoot->label = 1;
		-- labelCount[0];
		++ labelCount[1];


		++ numRelabels;


		addToStrongBucket (strongRoot, strongRoots[strongRoot->label].end);
	}	

	highestStrongLabel = 1;

	strongRoot = strongRoots[1].start->next;
	strongRoot->next->prev = strongRoot->prev;
	strongRoot->prev->next = strongRoot->next;
	strongRoot->next = NULL;

	return strongRoot;	
}



template <typename _Cap>
int optnet_pseudoflow<_Cap>::computeMinCut (void)
{
	size_type i, mincut = 0;

	for (i=0; i<numArc1s; ++i) 
	{
		/*
		if (Arc1List[i].from->label >=numNodes)
			labelList[Arc1List[i].from->number-1]=1;
		else
			labelList[Arc1List[i].from->number-1]=2;

		*/
		if ((Arc1List[i]->from->label >= numNodes) && (Arc1List[i]->to->label < numNodes))
		{
			mincut += Arc1List[i]->capacity;
		}
	}
	for (i=0;i<numNodes;i++)
	{
		if (adjacencyList[i].label>=numNodes)
			labelList[adjacencyList[i].number-1]=1;
		else
			labelList[adjacencyList[i].number-1]=2;

	}
	return mincut;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::pseudoflowPhase1 (void) 
{
	Node *strongRoot;
	int theparam = 0;
	

	while ((strongRoot = getHighestStrongRoot (theparam)))  
	{ 
		processRoot (strongRoot);
	}
	m_flow=computeMinCut();
	printf ("c Finished solving parameter %d\nc Flow: %ld\nc \n", 
		(theparam+1),
		m_flow);
	//std::cout<<"m_flow: "<<m_flow<<endl;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::checkOptimality (void) 
{
	int i, check = 1;
	llint mincut = 0, *excess; 

	//excess = (llint *) malloc (numNodes * sizeof (llint));
	excess= new llint [numNodes];
	/*
	if (!excess)
	{
		printf ("%s Line %d: Out of memory\n", __FILE__, __LINE__);
		exit (1);
	}
	*/

	for (i=0; i<numNodes; ++i)
	{
		excess[i] = 0;
	}

	for (i=0; i<numArc1s; ++i) 
	{
		if ((Arc1List[i]->from->label >= numNodes) && (Arc1List[i]->to->label < numNodes))
		{
			mincut += Arc1List[i]->capacity;
		}

		if ((Arc1List[i]->flow > Arc1List[i]->capacity) || (Arc1List[i]->flow < 0)) 
		{
			check = 0;
			printf("c Capacity constraint violated on Arc1 (%d, %d)\n", 
				Arc1List[i]->from->number,
				Arc1List[i]->to->number);
		}
		excess[Arc1List[i]->from->number - 1] -= Arc1List[i]->flow;
		excess[Arc1List[i]->to->number - 1] += Arc1List[i]->flow;
	}

	for (i=0; i<numNodes; i++) 
	{
		if ((i != (source-1)) && (i != (sink-1))) 
		{
			if (excess[i]) 
			{
				check = 0;
				printf ("c Flow balance constraint violated in node %d. Excess = %lld\n", 
					i+1,
					excess[i]);
			}
		}
	}

	if (check)
	{
		printf ("c\nc Solution checks as feasible.\n");
	}

	check = 1;

	if (excess[sink-1] != mincut) 
	{
		check = 0;
		printf("c Flow is not optimal - max flow does not equal min cut!\nc\n");
	}

	if (check) 
	{
		printf ("c\nc Solution checks as optimal.\nc \n");
		printf ("s Max Flow            : %lld\n", mincut);
	}

	free (excess);
	excess = NULL;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::quickSort (Arc1 **arr, const int first, const int last)
{
	int i, j, left=first, right=last, x1, x2, x3, mid, pivot, pivotval;
	Arc1 *swap;

	if ((right-left) <= 5)
	{// Bubble sort if 5 elements or less
		for (i=right; (i>left); --i)
		{
			swap = NULL;
			for (j=left; j<i; ++j)
			{
				if (arr[j]->flow < arr[j+1]->flow)
				{
					swap = arr[j];
					arr[j] = arr[j+1];
					arr[j+1] = swap;
				}
			}

			if (!swap)
			{
				return;
			}
		}

		return;
	}

	mid = (first+last)/2;

	x1 = arr[first]->flow; 
	x2 = arr[mid]->flow; 
	x3 = arr[last]->flow;

	pivot = mid;
	
	if (x1 <= x2)
	{
		if (x2 > x3)
		{
			pivot = left;

			if (x1 <= x3)
			{
				pivot = right;
			}
		}
	}
	else
	{
		if (x2 <= x3)
		{
			pivot = right;

			if (x1 <= x3)
			{
				pivot = left;
			}
		}
	}

	pivotval = arr[pivot]->flow;

	swap = arr[first];
	arr[first] = arr[pivot];
	arr[pivot] = swap;

	left = (first+1);

	while (left < right)
	{
		if (arr[left]->flow < pivotval)
		{
			swap = arr[left];
			arr[left] = arr[right];
			arr[right] = swap;
			-- right;
		}
		else
		{
			++ left;
		}
	}

	swap = arr[first];
	arr[first] = arr[left];
	arr[left] = swap;

	if (first < (left-1))
	{
		quickSort (arr, first, (left-1));
	}
	
	if ((left+1) < last)
	{
		quickSort (arr, (left+1), last);
	}
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::sort (Node * current)
{
	if (current->numOutOfTree > 1)
	{
		quickSort (current->outOfTree, 0, (current->numOutOfTree-1));
	}
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::minisort (Node *current) 
{
	Arc1 *temp = current->outOfTree[current->nextArc1];
	int i, size = current->numOutOfTree, tempflow = temp->flow;

	for(i=current->nextArc1+1; ((i<size) && (tempflow < current->outOfTree[i]->flow)); ++i)
	{
		current->outOfTree[i-1] = current->outOfTree[i];
	}
	current->outOfTree[i-1] = temp;
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::decompose (Node *excessNode, const int source, int *iteration) 
{
	Node *current = excessNode;
	Arc1 *tempArc1;
	int bottleneck = excessNode->excess;

	for ( ;(current->number != source) && (current->visited < (*iteration)); 
				current = tempArc1->from)
	{
		current->visited = (*iteration);
		tempArc1 = current->outOfTree[current->nextArc1];

		if (tempArc1->flow < bottleneck)
		{
			bottleneck = tempArc1->flow;
		}
	}

	if (current->number == source) 
	{
		excessNode->excess -= bottleneck;
		current = excessNode;

		while (current->number != source) 
		{
			tempArc1 = current->outOfTree[current->nextArc1];
			tempArc1->flow -= bottleneck;

			if (tempArc1->flow) 
			{
				minisort(current);
			}
			else 
			{
				++ current->nextArc1;
			}
			current = tempArc1->from;
		}
		return;
	}

	++ (*iteration);

	bottleneck = current->outOfTree[current->nextArc1]->flow;

	while (current->visited < (*iteration))
	{
		current->visited = (*iteration);
		tempArc1 = current->outOfTree[current->nextArc1];

		if (tempArc1->flow < bottleneck)
		{
			bottleneck = tempArc1->flow;
		}
		current = tempArc1->from;
	}	
	
	++ (*iteration);

	while (current->visited < (*iteration))
	{
		current->visited = (*iteration);

		tempArc1 = current->outOfTree[current->nextArc1];
		tempArc1->flow -= bottleneck;

		if (tempArc1->flow) 
		{
			minisort(current);
			current = tempArc1->from;
		}
		else 
		{
			++ current->nextArc1;
			current = tempArc1->from;
		}
	}
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::recoverFlow (void)
{
	int i, j, iteration = 1;
	Arc1 *tempArc1;
	Node *tempNode;

	for (i=0; i<adjacencyList[sink-1].numOutOfTree; ++i) 
	{
		tempArc1 = adjacencyList[sink-1].outOfTree[i];
		if (tempArc1->from->excess < 0) 
		{
			tempArc1->flow -= (int) (-1*tempArc1->from->excess); 
			tempArc1->from->excess = 0;
		}	
	}

	for (i=0; i<adjacencyList[source-1].numOutOfTree; ++i) 
	{
		tempArc1 = adjacencyList[source-1].outOfTree[i];
		addOutOfTreeNode (tempArc1->to, tempArc1);
	}

	adjacencyList[source-1].excess = 0;
	adjacencyList[sink-1].excess = 0;

	for (i=0; i<numNodes; ++i) 
	{
		tempNode = &adjacencyList[i];

		if ((i == (source-1)) || (i == (sink-1)))
		{
			continue;
		}

		if (tempNode->label >= numNodes) 
		{
			tempNode->nextArc1 = 0;
			if ((tempNode->parent) && (tempNode->Arc1ToParent->flow))
			{
				addOutOfTreeNode (tempNode->Arc1ToParent->to, tempNode->Arc1ToParent);
			}

			for (j=0; j<tempNode->numOutOfTree; ++j) 
			{
				if (!tempNode->outOfTree[j]->flow) 
				{
					-- tempNode->numOutOfTree;
					tempNode->outOfTree[j] = tempNode->outOfTree[tempNode->numOutOfTree];
					-- j;
				}
			}

			sort(tempNode);
		}
	}

	for (i=0; i<numNodes; ++i) 
	{
		tempNode = &adjacencyList[i];
		while (tempNode->excess > 0) 
		{
			++ iteration;
			decompose(tempNode, source, &iteration);
		}
	}
}

template <typename _Cap>
 void optnet_pseudoflow<_Cap>::displayBreakpoints (void)
{
	int i;
	for (i=0; i<numNodes; ++i)
	{
		printf ("n %d %d\n", (i+1), adjacencyList[i].breakpoint);
	}
}

template <typename _Cap>
void optnet_pseudoflow<_Cap>::freeMemory (void)
{
	int i;

	for (i=0; i<numNodes; ++i)
	{
		freeRoot (&strongRoots[i]);
	}

	free (strongRoots);

	for (i=0; i<numNodes; ++i)
	{
		if (adjacencyList[i].outOfTree)
		{
			free (adjacencyList[i].outOfTree);
		}
	}
	free (adjacencyList);
	free (labelCount);

}

} // namespace

#endif
