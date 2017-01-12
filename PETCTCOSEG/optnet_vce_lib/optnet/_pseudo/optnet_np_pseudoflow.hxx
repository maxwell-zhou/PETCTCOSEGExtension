/*
 ==========================================================================
 |   
 |   $Id: optnet_np_pseudoflow.hxx 21 2008-01-12 15:52:31Z Qi Song $
 |
 |   Written by Qi Song <qi-song@uiowa.edu>
 |   Department of Electrical and Computer Engineering
 |   University of Iowa
 |   
 ==========================================================================

 */

#ifndef ___OPTNET_PSEUDOFLOW_HXX___
#   define ___OPTNET_PSEUDOFLOW_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4127)
#   endif

//#   include <optnet/_alpha/graph_mc.hxx>
#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <queue>
////////////////////////////////////////////////////////////
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>


namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class optnet_np_pseudoflow_maxflow
///  @brief Implementation of the pseudoflow algorithm.
///////////////////////////////////////////////////////////////////////////
template <typename _Cap>
class optnet_pseudoflow
{
   

public:
 
	typedef _Cap  capacity_type;
    typedef unsigned int   size_type;
	

    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    optnet_pseudoflow();

    
    ///////////////////////////////////////////////////////////////////////
    ///  Construct a optnet_np_pseudoflow_maxflow object with the underlying graph
    ///  being created according to the given size information.
    ///
    ///  @param  size_z  The number of nodes in a column.
    ///  @param  size_x,size_y  Sizes for base graph.
    ///
    ///////////////////////////////////////////////////////////////////////
    bool create(size_type size_x, size_type size_y, size_type size_z, size_type size_s=1);

	bool create(size_type numpc, size_type numcols);

	///////////////////////////////////////////////////////////////////////
    ///  Return size information of x,y,z 
    ///////////////////////////////////////////////////////////////////////
	size_type size_0()
	{
		return m_x;
	}

	size_type size_1()
	{
		return m_y; 
	}

	size_type size_2()
	{
		return m_z;

	}
	size_type size_3()
	{
		return m_s;
	}
	size_type nodes_per_column()  //compatible with Kang's library, the same as m_z;
	{
   	return m_colsize;
	}
	size_type num_columns()  //compatible with Kang's library, the same as m_numcols;
	{
   	return m_numcols;
	}

    ///////////////////////////////////////////////////////////////////////
    ///  Solve the maximum-flow/minimum s-t cut problem.
    ///
    ///  @returns The maximum flow value.
    ///////////////////////////////////////////////////////////////////////
    capacity_type solve();
	

    ///////////////////////////////////////////////////////////////////////
    ///  Add arc(s) connecting a node to the source and/or the sink node.
    ///
    ///  @param  s    The capacity of the arc from the source node.
    ///  @param  t    The capacity of the arc to the sink node.
    ///  @param  index_x   The index of x
    ///  @param  index_y   The index of y
	///  @param  index_z   The index of z
    ///
    ///////////////////////////////////////////////////////////////////////
   void add_st_arc(capacity_type s, capacity_type t, size_type index_x, size_type index_y, size_type index_z, size_type index_s=0);

   void add_st_arc(capacity_type s, capacity_type t, size_type index_npc, size_type index_col);

   ///////////////////////////////////////////////////////////////////////
    ///  Add arc(s) connecting two nodes
    ///
    ///  @param  tail_x       start node
    ///  @param  tail_y       
    ///  @param  tail_z  
    ///  @param  head_x       end node
    ///  @param  head_y       
    ///  @param  head_z  
    ///
    ///////////////////////////////////////////////////////////////////////

   void add_arc(size_type tail_x, size_type tail_y, size_type tail_z, size_type head_x, size_type head_y, size_type head_z);

    ///////////////////////////////////////////////////////////////////////
    ///  Add arc(s) connecting two nodes
    ///
    ///  @param  tail_x       start node
    ///  @param  tail_y       
    ///  @param  tail_z 
    ///  @param  tail_s
    ///  @param  head_x       end node
    ///  @param  head_y       
    ///  @param  head_z  
    ///  @param  head_s
    ///////////////////////////////////////////////////////////////////////

   void add_arc(size_type tail_x, size_type tail_y, size_type tail_z, size_type tail_s, size_type head_x, size_type head_y, size_type head_z, size_type head_s);
   
   void add_arc(size_type tail_index_npc, size_type tail_index_col, size_type head_index_npc, size_type head_index_col);
    ///////////////////////////////////////////////////////////////////////
    ///  Add arc(s) connecting two nodes with edge cost
    ///  
    ///////////////////////////////////////////////////////////////////////
   void add_arc_cost(capacity_type edge_cost,size_type tail_index_npc, size_type tail_index_col, size_type head_index_npc, size_type head_index_col);
   
   void add_arc_cost(capacity_type edge_cost,size_type tail_x, size_type tail_y, size_type tail_z, size_type head_x, size_type head_y, size_type head_z);

   void add_arc_cost(capacity_type edge_cost,size_type tail_x, size_type tail_y, size_type tail_z, size_type tail_s, size_type head_x, size_type head_y, size_type head_z, size_type head_s);


    ///////////////////////////////////////////////////////////////////////
    ///  Determines if the given node is in the source set of the cut.
    ///
    ///  @param  index_x   
    ///  @param  index_y   
    ///  @param  index_z
    ///
    ///  @return Returns true if the given node is the source set,
    ///          false otherwise.
    ///
    ///////////////////////////////////////////////////////////////////////
    bool in_source_set(size_type index_x, size_type index_y, size_type index_z, size_type index_s=0);

	bool in_source_set(size_type index_npc, size_type index_col);

	/////////////////////////////////////////////////////////////////////////////
	void initialize(void);
	void freeMemory (void);
	////////////////////////////For compatibility, no use.
	void set_initial_flow(int tmpflow){};
	void clear_arcs(void){};



private:
	size_type m_colsize;
	size_type m_numcols;
	capacity_type m_flow;

	size_type m_x;
	size_type m_y;
	size_type m_z;
	size_type m_s;

	float m_lambda;

	typedef long long int llint;

    struct node;
    struct Arc1;

    typedef struct node 
     {
	  size_type visited;
	  size_type numAdjacent;
	  size_type number;
	  size_type label;
	  capacity_type excess;
	  struct node *parent;
	  struct node *childList;
	  struct node *nextScan;
	  size_type numOutOfTree;
	  Arc1 **outOfTree;
	  size_type nextArc1;
	  Arc1 *Arc1ToParent;
	  struct node *next;
	  struct node *prev;
	  int breakpoint;
     } Node;
	
	 typedef struct Arc1 
    {
	  struct node *from;
	  struct node *to;
	  capacity_type flow;
	  capacity_type capacity;
	  capacity_type direction;
	  capacity_type *capacities;
	  bool ispara;
    } Arc1;

    typedef struct root 
   {
	 Node *start;
	 Node *end;
   } Root;


    size_type numNodes;
    long numArc1s;
	long arcIndex;
    size_type source;
    size_type sink;
    size_type numParams;

    size_type highestStrongLabel;

	size_type *labelCount;
    Node *adjacencyList;
    Root *strongRoots;

	size_type *labelList;
    
   std::vector<Arc1*>  Arc1List;
	//std::deque<Arc1>  Arc1List;


  
     llint numPushes;
     int numMergers;
     int numRelabels;
     int numGaps;
     llint numArc1Scans;
    
     void prepareList();
	 void initializeNode (Node *nd, size_type n);
	 void initializeRoot (Root *rt);
	 void freeRoot (Root *rt);
	 void liftAll (Node *rootNode, size_type theparam);
	 void addToStrongBucket (Node *newRoot, Node *rootEnd);
	 void createOutOfTree (Node *nd);
	 void initializeArc1 (Arc1* ac);
	 void addOutOfTreeNode (Node *n, Arc1 *out);
	 void simpleInitialization (void);
	 int addRelationship (Node *newParent, Node *child);
	 void breakRelationship (Node *oldParent, Node *child);
	 void merge (Node *parent, Node *child, Arc1 *newArc1);
	 void pushUpward (Arc1 *currentArc1, Node *child, Node *parent, const int resCap);
	 void pushDownward (Arc1 *currentArc1, Node *child, Node *parent, int flow);
	 void pushExcess (Node *strongRoot);
	 Arc1 *findWeakNode (Node *strongNode, Node **weakNode);
	 void checkChildren (Node *curNode);
	 void processRoot (Node *strongRoot);
	 Node *getHighestStrongRoot (const int theparam);
	 int computeMinCut (void);
	 void pseudoflowPhase1 (void);
	 void checkOptimality (void);
	 void quickSort (Arc1 **arr, const int first, const int last);
	 void sort (Node * current);
	 void minisort (Node *current);
	 void decompose (Node *excessNode, const int source, int *iteration);
	 void recoverFlow (void);
	 void displayBreakpoints (void);

};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_pseudo/optnet_np_pseudoflow.cxx>
#   endif

#endif
