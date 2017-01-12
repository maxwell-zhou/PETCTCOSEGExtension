/*
 ==========================================================================
 |   Written by Qi Song <qi-song@uiowa.edu>
 |   Department of Electrical and Computer Engineering
 |   University of Iowa
 |   
 ==========================================================================
 */

#ifndef ___OPTNET_GS_GT_MULTI_DIR_HXX___
#   define ___OPTNET_GS_GT_MULTI_DIR_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_pseudo/optnet_np_pseudoflow.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <vector>


namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class optnet_gs_gt_multi_dir
///  @brief Implementation of the multiple-surface Optimal Net algorithm
///         using the Boykov-Kolmogorov max-flow algorithm on a
///         forward-star represented graph.
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg = net_f_xy>
class optnet_gs_gt_multi_dir
{
    //typedef optnet_fs_maxflow<_Cap, net_f_xy>   graph_type;
	typedef optnet_pseudoflow<_Cap>   graph_type;
/*
    struct  _Intra {
        std::vector<                    //
            std::pair<size_t, int>      // A tiny adjacency graph
        >       tonode[2];              //
        bool    circle[2];
        int     margin[2];
        int     smooth[2];
        int     visits[2];
        bool    dropped;
        bool    climbed;
    };
	*/

    struct  _Inter {
        size_t  k[2];
        int     r[2];
    };

	struct _Inter_cutsearch {      //Added by Sq
		size_t k[2];
		int r;
	};
	
	

	typedef std::vector<_Inter_cutsearch> inter_cutsearch_vector; //Added by Sq
    //typedef std::vector<_Intra>                 intra_vector;
    typedef std::vector<_Inter>                 inter_vector;
	
	struct _Shape_VCE{              //Added by Sq
		size_t	k;
		int		dir;				//Six directions, 0: x-positive; 1: x-negative;
                                    //                2: y-positive; 3: y-negative;
                                    //                4: z-positive; 5: z-negative; 									 
		int     **meanDir0;        
		int     **meanDir1;
		int     **upDir0;
		int     **upDir1;
		int     **lowDir0;
		int     **lowDir1;
		float	**fwdCofDir0;
		float	**backCofDir0;
		float	**fwdCofDir1;
		float	**backCofDir1;
	};

	typedef std::vector<_Shape_VCE>              shape_vce_vector;

public:
    typedef size_t                              size_type;
    typedef _Cost                               cost_type;
    typedef _Cap                                capacity_type;

    typedef array_base<cost_type, _Tg>          cost_array_base_type;
    typedef array_ref<cost_type, _Tg>           cost_array_ref_type;
    typedef array<cost_type, _Tg>               cost_array_type;
    
    typedef array_base<size_type>               net_base_type;
    typedef array_ref<size_type>                net_ref_type;
    typedef array<size_type>                    net_type; 
    typedef _Shape_VCE            				shape_vce_type;	
	//
	struct _Inter_cutcut {      //Added by Sq
		size_t k[2];
		cost_array_type* cost_context_cut;
	};
	//
	typedef _Inter_cutcut						inter_cutcut_type;
	typedef std::vector<_Inter_cutcut>          inter_cutcut_vector;
	//
	long flow_value;
	bool is_vce;
	int pow_vce;
	    
    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    optnet_gs_gt_multi_dir();

    ///////////////////////////////////////////////////////////////////////
    ///  Create graph and set the cost values.
    ///
    ///  @param cost The array that contain the cost values for each voxel.
    ///
    ///  @remarks This function will re-assign the node costs of the
    ///           underlying graph based upon the given cost array.
    ///
    ///////////////////////////////////////////////////////////////////////
    void create(const cost_array_base_type& cost);

	///////////////////////////////////////////////////////////////////////
    ///  Create graph for both graph search and graphcut subgraphs .
    ///
    ///  @param cost The array that contain the cost values for each voxel.
    ///
    ///  @remarks This function will re-assign the node costs of the
    ///           underlying graph based upon the given cost array.
    ///
    ///////////////////////////////////////////////////////////////////////
    void create( size_type s_0, size_type s_1, size_type s_2, size_type num_surf_graphsearch, size_type num_surf_graphcut );

    ///////////////////////////////////////////////////////////////////////
    ///  Set the smoothness constraints for the k-th surface in the graph 
	///  search framework.
    ///
    ///  @param k       The index to the k-th surface.
    ///  @param smooth0 The first smoothness parameter.
    ///  @param smooth1 The second smoothness parameter.
    ///  @param circle0 Enabling/disabling circle graph construction.
    ///  @param circle1 Enabling/disabling circle graph construction.
    ///
    ///  @remarks The actual meanings of the smoothness parameters depend
    ///           on the direction setting of the net. 
    ///
    ///////////////////////////////////////////////////////////////////////
    void set_graphsearch_params(size_type   k,
                    int         smooth0 = 1,
                    int         smooth1 = 1,
                    bool        circle0 = false,
                    bool        circle1 = false
                    );


    ///////////////////////////////////////////////////////////////////////
    ///  Set the relation between the surface in the graph cut
	///  framework and the surface in the graph search framework.
    ///
    ///  @param  k0  The index to one surface.
    ///  @param  k1  The index to the other surface.
    ///  @param  r01 The distance between the two surfaces.
    /// 
    ///
    ///////////////////////////////////////////////////////////////////////
    void set_cutsearch_relation(size_type k0, size_type k1, int r );
	
	///////////////////////////////////////////////////////////////////////
    ///  Set the relation between surfaces in the graph cut
	///  framework.
    ///
    ///  @param  k0  The index to one surface.
    ///  @param  k1  The index to the other surface.
    ///  @param  r01 The distance between the two surfaces.
    /// 
    ///
    ///////////////////////////////////////////////////////////////////////
    void set_cutcut_relation(const inter_cutcut_type& context_cut ){ m_inter_cutcut.push_back( context_cut );};
	
	///////////////////////////////////////////////////////////////////////
    ///  Set the shape prior.
    ///
    ///  @param  shape_vce
    ///
    ///////////////////////////////////////////////////////////////////////
	void set_shape_prior(const shape_vce_type& shape_vce);

    ///////////////////////////////////////////////////////////////////////
    ///  Solve the optimal surface problem using the given cost function
    ///  and smoothness constraints.
    ///
    ///  @param net   The resulting optimal "net" surface.
    ///  @param pflow The output maximum flow value.
    ///
    ///  @remarks The pflow parameter, if not NULL, will return the
    ///           computed maximum-flow value. It is used primarily
    ///           for debugging.
    ///
    ///////////////////////////////////////////////////////////////////////
    void solve(net_base_type& net,      // [OUT]
               capacity_type* pflow = 0 // [OUT]
               );

	///////////////////////////////////////////////////////////////////////
    ///  Find the optimal cut.
    ///
    ///  @param net   The resulting labeled image..
    ///  @param pflow The output maximum flow value.
    ///
    ///  @remarks The pflow parameter, if not NULL, will return the
    ///           computed maximum-flow value. It is used primarily
    ///           for debugging.
    ///
    ///////////////////////////////////////////////////////////////////////
    void solve_all (net_base_type& net,      // [OUT]
               capacity_type* pflow = 0 // [OUT]
               );
	///////////////////////////////////////////////////////////////////////
	// Set cost of nodes in graph search framework.
	void set_gs_cost(const cost_array_type& cost) { m_pcost_gs = &cost; };

	///////////////////////////////////////////////////////////////////////
	// Set cost of nodes belong to object (regional term of energy of graph cut)
	void set_ob_cost(const cost_array_type& cost) { m_pcost_ob = &cost; };

	///////////////////////////////////////////////////////////////////////
	// Set cost of nodes belong to background (regional term of energy of graph cut)
	void set_bg_cost(const cost_array_type& cost) { m_pcost_bg = &cost; };

	///////////////////////////////////////////////////////////////////////
	// Set cost of pairs of neighboring nodes (boundary term of energy of graph cut)
	void set_neigh_cost(const cost_array_type& cost) { m_pcost_neigh = &cost; };
	
	void set_neigh_coef(capacity_type* coef) { m_neigh_coef = coef; }
	
	
private:

    ///////////////////////////////////////////////////////////////////////
    // Compute the upper and lower margin of the 3-D subgraphs. The nodes
    // above the upper bound and below the lower bound can never be on
    // the resulting surfaces.
    //void get_bounds_of_subgraphs();

    ///////////////////////////////////////////////////////////////////////
    // Transform the costs of the graph nodes based on the given
    // cost vector.
    void transform_costs_x(const shape_vce_type& shape_vce);
	
	void transform_costs_y(const shape_vce_type& shape_vce);
	
	void transform_costs_z(const shape_vce_type& shape_vce);
    
    ///////////////////////////////////////////////////////////////////////
    // Construct the arcs of the underlying graph.
    //void build_arcs();

	///////////////////////////////////////////////////////////////////////
    // Construct the vce arcs of the underlying graph in three directions.
    void build_vce_arcs_x(const shape_vce_type& shape_vce);
	
	void build_vce_arcs_y(const shape_vce_type& shape_vce);
	
	void build_vce_arcs_z(const shape_vce_type& shape_vce);
	

	///////////////////////////////////////////////////////////////////////
    // Construct the graph using grachcut methd.
    void build_graphcut_arcs();

	////////////////////////////////////////////////////////////////////////
	// Construct inter-surface arcs between sub-graphs of graph cut and graph 
	// search
    void build_gs_gc_arcs();
	
	////////////////////////////////////////////////////////////////////////
	// Construct inter-surface arcs between sub-graphs of graph cut
    void build_gc_gc_arcs();

    ///////////////////////////////////////////////////////////////////////
	// Pointer for cost of nodes (graph search)
    const cost_array_base_type* m_pcost_gs;

	///////////////////////////////////////////////////////////////////////
	// Pointer for cost of nodes belonging to object (regional term of 
	// energy of graph cut).
    const cost_array_type* m_pcost_ob;

	///////////////////////////////////////////////////////////////////////
	// Pointer for cost of nodes belonging to background (regional term of 
	// energy of graph cut).
    const cost_array_type* m_pcost_bg;

	///////////////////////////////////////////////////////////////////////
	// Pointer for cost of pairs of neighboring nodes (boundry term of 
	// energy of graph cut)
    const cost_array_type* m_pcost_neigh;
	
	///////////////////////////////////////////////////////////////////////
	//int arc_weight( int k , size_type i0, size_type i1, size_type i3, int dir, int fwdFlag);
	
	///////////////////////////////////////////////////////////////////////
	int arc_weight( int k );

    graph_type                m_graph;
    //intra_vector              m_intra;
    inter_vector              m_inter;
	inter_cutsearch_vector    m_inter_cutsearch;
	inter_cutcut_vector		  m_inter_cutcut;
	shape_vce_vector          m_shape_prior;
	size_type                 m_num_surf_graphsearch;
	size_type                 m_num_surf_graphcut;
	capacity_type*			  m_neigh_coef;
};

} // optnet

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet_graphcut/optnet_gs_gt_multi_dir.cxx>
#   endif

#endif
