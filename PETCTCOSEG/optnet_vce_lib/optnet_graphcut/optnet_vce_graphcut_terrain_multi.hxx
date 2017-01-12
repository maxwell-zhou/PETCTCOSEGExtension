/*
 ==========================================================================
 |   Written by Qi Song <qi-song@uiowa.edu>
 |   Department of Electrical and Computer Engineering
 |   University of Iowa
 |   
 ==========================================================================
 */

#ifndef ___OPTNET_vce_3d_terrain_MULTI_HXX___
#   define ___OPTNET_vce_3d_terrain_MULTI_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
//#   include <optnet/_fs/optnet_fs_maxflow.hxx>
#   include <optnet/_pseudo/optnet_np_pseudoflow.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <vector>


namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class optnet_vce_graphcut_terrain_multi
///  @brief Implementation of the multiple-surface Optimal Net algorithm
///         using the Boykov-Kolmogorov max-flow algorithm on a
///         forward-star represented graph.
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg = net_f_xy>
class optnet_vce_graphcut_terrain_multi
{
    //typedef optnet_fs_maxflow<_Cap, net_f_xy>   graph_type;
	typedef optnet_pseudoflow<_Cap>   graph_type;

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

    struct  _Inter {
        size_t  k[2];
        int     r[2];
    };

	struct _Inter_cutsearch {      //Added by Sq
		size_t k[2];
		int r;
	};

	typedef std::vector<_Inter_cutsearch> inter_cutsearch_vector; //Added by Sq

    typedef std::vector<_Intra>                 intra_vector;
    typedef std::vector<_Inter>                 inter_vector;
	
	struct _Smooth{               //Added by Sq
		int     **smoothDir0;        
		int     **smoothDir1;        
	};

	typedef std::vector<_Smooth>           smooth_vector;

	struct _ArcCof{               //Added by Sq
		int     **fwdDir0;
		int     **backDir0;
		int     **fwdDir1;
		int     **backDir1;
	};

	typedef std::vector<_ArcCof>              arc_cof_vector;


	struct _VCE{               //Added by Sq
		int     **meanDir0;        
		int     **meanDir1;
		int     **upDir0;
		int     **upDir1;
		int     **lowDir0;
		int     **lowDir1;
	};

	typedef std::vector<_VCE>              vce_vector;

	


public:

	/*
    typedef typename graph_type::roi_base_type  roi_base_type;
    typedef typename graph_type::roi_ref_type   roi_ref_type;
    typedef typename graph_type::roi_type       roi_type;
	*/

    typedef size_t                              size_type;
    typedef _Cost                               cost_type;
    typedef _Cap                                capacity_type;

    typedef array_base<cost_type, _Tg>          cost_array_base_type;
    typedef array_ref<cost_type, _Tg>           cost_array_ref_type;
    typedef array<cost_type, _Tg>               cost_array_type;
    
    typedef array_base<size_type>               net_base_type;
    typedef array_ref<size_type>                net_ref_type;
    typedef array<size_type>                    net_type;
    
	//
	
	long flow_value;
	bool is_vce;

	capacity_type **arc_cost_pos;
	capacity_type **arc_cost_neg;
	int pow_vce;
	smooth_vector   smooth_array;
	vce_vector      vce_para;
	arc_cof_vector  arc_cof;
    
    ///////////////////////////////////////////////////////////////////////
    /// Default constructor.
    ///////////////////////////////////////////////////////////////////////
    optnet_vce_graphcut_terrain_multi();

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


private:

    ///////////////////////////////////////////////////////////////////////
    // Compute the upper and lower margin of the 3-D subgraphs. The nodes
    // above the upper bound and below the lower bound can never be on
    // the resulting surfaces.
    void get_bounds_of_subgraphs();

    ///////////////////////////////////////////////////////////////////////
    // Transform the costs of the graph nodes based on the given
    // cost vector.
    void transform_costs();
    
    ///////////////////////////////////////////////////////////////////////
    // Construct the arcs of the underlying graph.
    void build_arcs();

	///////////////////////////////////////////////////////////////////////
    // Construct the vce arcs of the underlying graph.
    void build_vce_arcs();

	///////////////////////////////////////////////////////////////////////
    // Construct the graph using grachcut methd.
    void build_graphcut_arcs();

	////////////////////////////////////////////////////////////////////////
	// Construct inter-surface arcs
    void build_intersurface_arcs();

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
	int arc_weight( int k , size_type i0, size_type i1, size_type i3, int dir, int fwdFlag);

    graph_type                m_graph;
    intra_vector              m_intra;
    inter_vector              m_inter;
	inter_cutsearch_vector    m_inter_cutsearch;
	size_type                 m_num_surf_graphsearch;
	size_type                 m_num_surf_graphcut;
};

} // optnet

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet_graphcut/optnet_vce_graphcut_terrain_multi.cxx>
#   endif

#endif
