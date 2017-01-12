/*
 ==========================================================================
 |   Written by Qi Song <qi-song@uiowa.edu>
 |   Department of Electrical and Computer Engineering
 |   University of Iowa
 |   
 ==========================================================================
 */

#ifndef ___OPTNET_GS_GT_MULTI_DIR_CXX___
#   define ___OPTNET_GS_GT_MULTI_DIR_CXX___

#   include <optnet/_base/except.hxx>
#   include <optnet_graphcut/optnet_gs_gt_multi_dir.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <deque>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::optnet_gs_gt_multi_dir() :
    m_pcost_gs(0)
{}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::create( size_type s_0, size_type s_1, size_type s_2, size_type num_surf_graphsearch, size_type num_surf_graphcut )
{
	size_type s_3 = num_surf_graphsearch + num_surf_graphcut;

    if ( !m_graph.create( s_0, s_1, s_2, s_3 ) ) 
	{
        throw_exception(std::runtime_error(
            "optnet_gs_gt_multi_dir::create: Could not create graph."
        ));
    }

    m_num_surf_graphsearch = num_surf_graphsearch;
	m_num_surf_graphcut = num_surf_graphcut;

    m_inter.clear();

	pow_vce = 1;

}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
int
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::arc_weight(int k)
                                          
{
   int weight;
   if ( k == 0 )
   weight = 1;
   else
   weight = pow( double(k+1), double(pow_vce) )-2 * pow( double(k), double(pow_vce)) + pow( double(k-1), double(pow_vce)) ;
   
   return weight;
   
}


//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::set_cutsearch_relation(size_type k0,
                                                   size_type k1,
                                                   int       r
                                                   )
{
 
    // Create a new relation and 
    //  append it to the relation vector.
    _Inter_cutsearch relation;

    relation.r = r;
    relation.k[0] = k0;
	relation.k[1] = k1;
    
    m_inter_cutsearch.push_back( relation );

}
//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::set_shape_prior(const shape_vce_type& shape_vce)
{    
    m_shape_prior.push_back( shape_vce );
}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::solve_all(net_base_type& net, 
                                            capacity_type* pflow
                                            )
{
    size_type       i0, i1, i2, i3;
	int 			s0, s1, s2;
    capacity_type   flow;

    if (net.size_0() != m_graph.size_0() || 
        net.size_1() != m_graph.size_1() ||
		net.size_2() != m_graph.size_2() ||
        net.size_3() != m_graph.size_3()
        ) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_gs_gt_multi_dir::solve: The output image size must match the graph size."
        ));
    }

	m_graph.set_initial_flow(0);
	// Assign the cost of graph nodes based on the input cost
    // vector. We also perform the "translation operation"
    // here to guaranttee a non-empty solution.
   
    //Build the graph for graph search.
    for (i3 = 0; i3 < m_shape_prior.size(); i3++)
	{
		if ( m_shape_prior[i3].dir == 0 || m_shape_prior[i3].dir == 1 )
		{
			transform_costs_x( m_shape_prior[i3] );
			build_vce_arcs_x( m_shape_prior[i3] );
		}
		else if ( m_shape_prior[i3].dir == 2 || m_shape_prior[i3].dir == 3 )
		{
			transform_costs_y( m_shape_prior[i3] );
			build_vce_arcs_y( m_shape_prior[i3] );
		}
		else if ( m_shape_prior[i3].dir == 4 || m_shape_prior[i3].dir == 5 )
		{
			transform_costs_z( m_shape_prior[i3] );
			build_vce_arcs_z( m_shape_prior[i3] );
		}
		
	}
	//transform_costs();
    //cout<<"Finish cost transform"<<endl;

    // Build the arcs of the graphs.
	//build_vce_arcs();
	cout << "Build graph cut arcs" << endl;
    build_graphcut_arcs();
	cout << "Build gs gc arcs" << endl;
	build_gs_gc_arcs();
	cout << "Build gc gc arcs" << endl;
	build_gc_gc_arcs();
	cout << "Finish build arcs" << endl;

    // Calculate max-flow/min-cut.
    flow = m_graph.solve();
	
	//Get the labeled image for graph search.
    for (i3 = 0; i3 < m_shape_prior.size(); i3++)
	{

		
		if ( m_shape_prior[i3].dir == 0 || m_shape_prior[i3].dir == 1 )
		{
			int graph_id = m_shape_prior[i3].k;
			s0 = (int)m_graph.size_1();
			s1 = (int)m_graph.size_2();
			s2 = (int)m_graph.size_0();
			for (i1 = 0; i1 < s1; ++i1)
                for (i0 = 0; i0 < s0; ++i0) 
				{
					for ( i2 = 0; i2 < s2; ++i2) 
					{
                        // Find upper envelope.
                        if (!m_graph.in_source_set(i2, i0, i1, graph_id))
                             break;
                    }
					if ( m_shape_prior[i3].dir == 0)
						net( i2, i0, i1, graph_id ) = 1;
					else 
						net( s2 - i2 - 1, i0, i1, graph_id ) = 1;		
				}
		}
		else if ( m_shape_prior[i3].dir == 2 || m_shape_prior[i3].dir == 3 )
		{
			int graph_id = m_shape_prior[i3].k;
			s0 = (int)m_graph.size_0();
			s1 = (int)m_graph.size_2();
			s2 = (int)m_graph.size_1();
			for (i1 = 0; i1 < s1; ++i1)
                for (i0 = 0; i0 < s0; ++i0) 
				{
					for ( i2 = 0; i2 < s2; ++i2) 
					{
                        // Find upper envelope.
                        if (!m_graph.in_source_set(i0, i2, i1, graph_id))
                             break;
                    }
					if ( m_shape_prior[i3].dir == 2)
						net( i0, i2, i1, graph_id ) = 1;
					else 
						net( i0, s2 - i2 - 1, i1, graph_id ) = 1;		
				}
		}
		else if ( m_shape_prior[i3].dir == 4 || m_shape_prior[i3].dir == 5 )
		{
			int graph_id = m_shape_prior[i3].k;
			s0 = (int)m_graph.size_0();
			s1 = (int)m_graph.size_1();
			s2 = (int)m_graph.size_2();
			for (i1 = 0; i1 < s1; ++i1)
                for (i0 = 0; i0 < s0; ++i0) 
				{
					for ( i2 = 0; i2 < s2; ++i2) 
					{
                        // Find upper envelope.
                        if (!m_graph.in_source_set(i0, i1, i2, graph_id))
                             break;
                    }
					if ( m_shape_prior[i3].dir == 4)
						net( i0, i1, i2, graph_id ) = 1;
					else 
						net( i0, i1, s2 - i2 - 1, graph_id ) = 1;		
				}
			std::cout << "Compute the resulted image" << std::endl;
		}
		
	}
	
	//Get the labeled image for graph cut.
	for ( i3 = m_num_surf_graphsearch; i3 < m_graph.size_3(); ++i3)
	{
		for (i1 = 0; i1 < m_graph.size_1(); ++i1) 
            for (i0 = 0; i0 < m_graph.size_0(); ++i0) 
				for (i2 = 0; i2 < m_graph.size_2(); ++i2)
				{
                    if ( m_graph.in_source_set( i0, i1, i2, i3 ) )
						net( i0, i1, i2, i3 ) = 1;
					else
						net( i0, i1, i2, i3 ) = 0;	
				}					
	}
    
    if (0 != pflow)
        *pflow = flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::transform_costs_z(const shape_vce_type& shape_vce)
{
    assert(m_pcost_gs != 0);
    size_type            i0, i1, i2, i3, s0, s1, s2, s3;

    s0 = m_pcost_gs->size_0();
    s1 = m_pcost_gs->size_1();
    s2 = m_pcost_gs->size_2();
    //s3 = m_pcost_gs->size_3();
	i3 = shape_vce.k;
    //m_graph.set_initial_flow(0);
    // Construct the s-t graph "G_st".

    for (i2 = 1;	i2 < s2;	++i2) 
	{
        for (i1 = 0; i1 < s1; ++i1) 
		{
            for (i0 = 0; i0 < s0; ++i0) 
			{
				capacity_type cap;
                   
				if ( shape_vce.dir == 4 )
                    cap = (capacity_type)(*m_pcost_gs)(i0, i1, i2, i3) - (capacity_type)(*m_pcost_gs)(i0, i1, i2 - 1, i3);
				else     //Reverse the image. 
					cap = (capacity_type)(*m_pcost_gs)(i0, i1, s2 - i2 - 1, i3) - (capacity_type)(*m_pcost_gs)(i0, i1, s2 - i2, i3);

                if (cap >= 0) // non-negative -> connect to t
                    m_graph.add_st_arc(0, +cap, i0, i1, i2, i3);
                else          // negative     -> connect to s
                    m_graph.add_st_arc(-cap, 0, i0, i1, i2, i3);

            } // for i0
        } // for i1
    } // for i2
 
	for (i1 = 0; i1 < s1; ++i1) 
	{
        for (i0 = 0; i0 < s0; ++i0)
		{                    
            m_graph.add_st_arc(100000, 0, i0, i1, 0, i3);
        } // for i0
    } // for i1
}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::transform_costs_x(const shape_vce_type& shape_vce)
{
    assert(m_pcost_gs != 0);
    size_type            i0, i1, i2, i3, s0, s1, s2, s3;

    s0 = m_pcost_gs->size_1();
    s1 = m_pcost_gs->size_2();
    s2 = m_pcost_gs->size_0();
    
	i3 = shape_vce.k;
    //m_graph.set_initial_flow(0);
    // Construct the s-t graph "G_st".

    for (i2 = 1;	i2 < s2;	++i2) 
	{
        for (i1 = 0; i1 < s1; ++i1) 
		{
            for (i0 = 0; i0 < s0; ++i0) 
			{
				capacity_type cap;
                   
				if ( shape_vce.dir == 0 )
                    cap = (capacity_type)(*m_pcost_gs)(i2, i0, i1, i3) - (capacity_type)(*m_pcost_gs)(i2 - 1, i0, i1, i3);
				else     //Reverse the image. 
					cap = (capacity_type)(*m_pcost_gs)(s2 - i2 - 1, i0, i1,  i3) - (capacity_type)(*m_pcost_gs)(s2 - i2, i0, i1,  i3);

                if (cap >= 0) // non-negative -> connect to t
                    m_graph.add_st_arc(0, +cap, i2, i0, i1, i3);
                else          // negative     -> connect to s
                    m_graph.add_st_arc(-cap, 0, i2, i0, i1, i3);

            } // for i0
        } // for i1
    } // for i2
 
	for (i1 = 0; i1 < s1; ++i1) 
	{
        for (i0 = 0; i0 < s0; ++i0)
		{                    
            m_graph.add_st_arc(100000, 0, 0, i0, i1, i3);
        } // for i0
    } // for i1
}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::transform_costs_y(const shape_vce_type& shape_vce)
{
    assert(m_pcost_gs != 0);
    size_type            i0, i1, i2, i3, s0, s1, s2, s3;

    s0 = m_pcost_gs->size_0();
    s1 = m_pcost_gs->size_2();
    s2 = m_pcost_gs->size_1();
    
	i3 = shape_vce.k;
    //m_graph.set_initial_flow(0);
    // Construct the s-t graph "G_st".

    for (i2 = 1;	i2 < s2;	++i2) 
	{
        for (i1 = 0; i1 < s1; ++i1) 
		{
            for (i0 = 0; i0 < s0; ++i0) 
			{
				capacity_type cap;
                   
				if ( shape_vce.dir == 2 )
                    cap = (capacity_type)(*m_pcost_gs)(i0, i2, i1, i3) - (capacity_type)(*m_pcost_gs)(i0, i2 - 1, i1, i3);
				else     //Reverse the image. 
					cap = (capacity_type)(*m_pcost_gs)(i0, s2 - i2 - 1, i1,  i3) - (capacity_type)(*m_pcost_gs)(i0, s2 - i2, i1,  i3);

                if (cap >= 0) // non-negative -> connect to t
                    m_graph.add_st_arc(0, +cap, i0, i2, i1, i3);
                else          // negative     -> connect to s
                    m_graph.add_st_arc(-cap, 0, i0, i2, i1, i3);

            } // for i0
        } // for i1
    } // for i2
 
	for (i1 = 0; i1 < s1; ++i1) 
	{
        for (i0 = 0; i0 < s0; ++i0)
		{                    
            m_graph.add_st_arc(100000, 0, i0, 0, i1, i3);
        } // for i0
    } // for i1
}

/////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::build_vce_arcs_z(const shape_vce_type& shape_vce)
{
	cout<<"Begin vce arcs building"<<endl;
    int i0, i1, i2, i3, s0, s1, s2, s3, ii;
    int convexPower;
	//int graph_id;
	
	i3 = shape_vce.k;
	convexPower = pow_vce;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_1();
    s2 = (int)m_graph.size_2();
	//s3 = (int)(m_num_surf_graphsearch);

    // Construct graph arcs based on the given parameters.
    const int & bounds0 = 0;
    const int & bounds1 = 0;

	//// -- Intra-column (vertical) arcs.
    for (i2 = s2 - bounds1 - 1;
             i2 > bounds0;
             --i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1,          i2 - 1,          i3);
            } // for i0 
        } // for i1
    } // for i2

    // Inter-column arcs (dir-0).

	if ( i0 > 1)
	{
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0 - 1; ++i0) {
				
				int upBound = 0, lowBound = 0; 
				//Forward
                if ( shape_vce.meanDir0[i0][i1] > 0 )
				upBound = shape_vce.meanDir0[i0][i1];
				//cout<<"The upBound is: "<<upBound<<endl;

				if ( (shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]) < 0 )
				lowBound = - ( shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]);

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.lowDir0[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.fwdCofDir0[i0][i1] * arc_weight(k),     i0,          i1,          i2,          i3,		i0 + 1,      i1,          i2 + shape_vce.meanDir0[i0][i1] - k,		i3 );
					
					m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1],		i3 );//Hard constraint
				
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
				{
					if ( (i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1],		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          s2 - 1,		i3 );//Hard constraint
				
					
				
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]) >= 0 )
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1],		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          0,		i3 );//Hard constraint
				}
				//End boundary condition
				//Backward
				upBound = 0, lowBound = 0;
				if ( shape_vce.meanDir0[i0][i1] < 0 )
				upBound = -shape_vce.meanDir0[i0][i1];

				if ( (shape_vce.meanDir0[i0][i1] + shape_vce.upDir0[i0][i1]) > 0 )
				lowBound = shape_vce.meanDir0[i0][i1] + shape_vce.upDir0[i0][i1];

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.upDir0[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.backCofDir0[i0][i1] * arc_weight(k),     i0 + 1,          i1,          i2,          i3,		i0,      i1,          i2 - shape_vce.meanDir0[i0][i1] - k,		i3 );

					m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1],		i3 );//Hard constraints
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
					{
					if ( (i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1],		i3 );//Hard constraints
					else
						m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          s2 - 1,		i3 );//Hard constraints
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1]) >= 0 )
						m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1],		i3 );//Hard constraints
					else
						m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          0,		i3 );//Hard constraints
				}
				//End boundary condition		
			} //for i0
		}// for i1
	} // for if
   cout<<"Finish dir-0"<<endl;
    // Inter-column arcs (dir-1). For 2-D case, no need
   if ( i1 > 1)
	{
    for (i1 = 0; i1 < s1 - 1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
				
				int upBound = 0, lowBound = 0; 
				//Forward
                if ( shape_vce.meanDir1[i0][i1] > 0 )
				upBound = shape_vce.meanDir1[i0][i1];

				if ( (shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]) < 0 )
				lowBound = - ( shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]);

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.lowDir1[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.fwdCofDir1[i0][i1] * arc_weight(k),     i0,          i1,          i2,          i3,		i0,      i1 + 1,          i2 + shape_vce.meanDir1[i0][i1] - k,		i3 );
					
					m_graph.add_arc( i0,          i1,          i2,          i3,		i0,      i1 + 1,          i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1],		i3 );//Hard constraint
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
				{
					if ( (i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0,      i1 + 1,          i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1],		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0,      i1 + 1,          s2 - 1,		i3 );//Hard constraint
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]) >= 0 )
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0,      i1 + 1,          i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1],		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0,      i1 + 1,          0,		i3 );//Hard constraint
				}
				//End boundary condition

				//Backward
				upBound = 0, lowBound = 0;
				if ( shape_vce.meanDir1[i0][i1] < 0 )
				upBound = -shape_vce.meanDir1[i0][i1];

				if ( (shape_vce.meanDir1[i0][i1] + shape_vce.upDir1[i0][i1]) > 0 )
				lowBound = shape_vce.meanDir1[i0][i1] + shape_vce.upDir1[i0][i1];

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.upDir1[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.backCofDir1[i0][i1] * arc_weight(k),     i0,          i1 + 1,          i2,          i3,		i0,      i1,          i2 - shape_vce.meanDir1[i0][i1] - k,		i3 );
				    
					m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1],		i3 ); //Hard constraint
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
					{
					if ( (i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1],		i3 );//Hard constraints
					else
						m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          s2 - 1,		i3 );//Hard constraints
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1]) >= 0 )
						m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1],		i3 );//Hard constraints
					else
						m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          0,		i3 );//Hard constraints
				}
				//End boundary condition		
			} //for i0
		}// for i1
	}// for if
    
    cout<<"Finish dir-1"<<endl;
    //Free memory
		//vce_para.clear();
		//arc_cof.clear();
   
} //build vce_z

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::build_vce_arcs_x(const shape_vce_type& shape_vce)
{
	cout<<"Begin vce arcs building"<<endl;
    int i0, i1, i2, i3, s0, s1, s2, s3, ii;
    int convexPower;
	int dir;
	
	i3 = shape_vce.k;
	convexPower = pow_vce;
	
	s0 = (int)m_graph.size_1();
	s1 = (int)m_graph.size_2();
	s2 = (int)m_graph.size_0();

	//s3 = (int)(m_num_surf_graphsearch);
    cout << "Build vce arcs x" << endl;
    // Construct graph arcs based on the given parameters.
    const int & bounds0 = 0;
    const int & bounds1 = 0;

	//// -- Intra-column (vertical) arcs.
    for (i2 = s2 - bounds1 - 1;
             i2 > bounds0;
             --i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
						m_graph.add_arc(i2,          i0,          i1,          i3,          i2 - 1,          i0,          i1,          i3);
            } // for i0 
        } // for i1
    } // for i2

    // Inter-column arcs (dir-0).
	if ( i0 > 1)
	{
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0 - 1; ++i0) {
				int upBound = 0, lowBound = 0; 
				//Forward
                if ( shape_vce.meanDir0[i0][i1] > 0 )
				upBound = shape_vce.meanDir0[i0][i1];
				
				if ( (shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]) < 0 )
				lowBound = - ( shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]);

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.lowDir0[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.fwdCofDir0[i0][i1] * arc_weight(k),     i2,          i0,          i1,          i3,		i2 + shape_vce.meanDir0[i0][i1] - k,      i0 + 1,          i2,		i3 );
						
					m_graph.add_arc(i2,          i0,          i1,          i3,		i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1],      i0 + 1,          i1,		i3 );//Hard constraint
				
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
				{
					if ( (i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc(i2,          i0,          i1,          i3,		i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1],      i0 + 1,          i1,		i3 );//Hard constraint
					else
						m_graph.add_arc(i2,          i0,          i1,          i3,		s2 - 1,      i0 + 1,     i1,		i3 );//Hard constraint
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]) >= 0 )
						m_graph.add_arc(i2,          i0,          i1,          i3,		i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1],      i0 + 1,          i1,		i3 );//Hard constraint
					else
						m_graph.add_arc(i2,          i0,          i1,          i3,		0,      i0 + 1,     i1,		i3 );//Hard constraint
				}
				//End boundary condition
				//Backward
				upBound = 0, lowBound = 0;
				if ( shape_vce.meanDir0[i0][i1] < 0 )
				upBound = -shape_vce.meanDir0[i0][i1];

				if ( (shape_vce.meanDir0[i0][i1] + shape_vce.upDir0[i0][i1]) > 0 )
				lowBound = shape_vce.meanDir0[i0][i1] + shape_vce.upDir0[i0][i1];
                
				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.upDir0[i0][i1]; k++)
				    {m_graph.add_arc_cost(shape_vce.backCofDir0[i0][i1] * arc_weight(k),     i2,          i0 + 1,          i1,          i3,		i2 - shape_vce.meanDir0[i0][i1] - k,      i0,          i1,		i3 );
					//cout << "Wrong detected!" << endl;
					}

					m_graph.add_arc( i2,          i0 + 1,          i1,          i3,		i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1],      i0,          i1,		i3 );//Hard constraints
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
					{
					if ( (i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc( i2,          i0 + 1,          i1,          i3,		i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1],      i0,          i1,		i3 );//Hard constraints
					else
						m_graph.add_arc( i2,          i0 + 1,          i1,          i3,		s2 - 1,      i0,          i1,		i3 );//Hard constraints
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1]) >= 0 )
						m_graph.add_arc( i2,          i0 + 1,          i1,          i3,		i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1],      i0,          i1,		i3 );//Hard constraints
					else
						m_graph.add_arc( i2,          i0 + 1,          i1,          i3,		0,      i0,     i1,		i3 );//Hard constraints
				}
				//End boundary condition		
			} //for i0
		}// for i1
	} // for if
   cout<<"Finish dir-0"<<endl;
    // Inter-column arcs (dir-1). For 2-D case, no need
   if ( i1 > 1)
	{
    for (i1 = 0; i1 < s1 - 1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
				
				int upBound = 0, lowBound = 0; 
				//Forward
                if ( shape_vce.meanDir1[i0][i1] > 0 )
				upBound = shape_vce.meanDir1[i0][i1];

				if ( (shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]) < 0 )
				lowBound = - ( shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]);

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.lowDir1[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.fwdCofDir1[i0][i1] * arc_weight(k),     i2,          i0,          i1,          i3,		i2 + shape_vce.meanDir1[i0][i1] - k,      i0,          i1 + 1,		i3 );
					
					m_graph.add_arc( i2,          i0,          i1,          i3,		i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1],      i0,          i1 + 1,		i3 );//Hard constraint
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
				{
					if ( (i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc(i2,          i0,          i1,          i3,		i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1],      i0,          i1 + 1,		i3 );//Hard constraint
					else
						m_graph.add_arc(i2,          i0,          i1,          i3,		s2 - 1,      i0,          i1 + 1,		i3 );//Hard constraint
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]) >= 0 )
						m_graph.add_arc(i2,          i0,          i1,          i3,		i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1],      i0,          i1 + 1,		i3 );//Hard constraint
					else
						m_graph.add_arc(i2,          i0,          i1,          i3,		0,      i0,          i1 + 1,		i3 );//Hard constraint
				}
				//End boundary condition

				//Backward
				upBound = 0, lowBound = 0;
				if ( shape_vce.meanDir1[i0][i1] < 0 )
				upBound = -shape_vce.meanDir1[i0][i1];

				if ( (shape_vce.meanDir1[i0][i1] + shape_vce.upDir1[i0][i1]) > 0 )
				lowBound = shape_vce.meanDir1[i0][i1] + shape_vce.upDir1[i0][i1];

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.upDir1[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.backCofDir1[i0][i1] * arc_weight(k),     i2,          i0,          i1 + 1,          i3,		i2 - shape_vce.meanDir1[i0][i1] - k,      i0,          i1,		i3 );
				    
					m_graph.add_arc( i2,          i0,          i1 + 1,          i3,		i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1],      i0,          i1,		i3); //Hard constraint
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
					{
					if ( (i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc( i2,          i0,          i1 + 1,          i3,		i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1],      i0,          i1,		i3 );//Hard constraints
					else
						m_graph.add_arc( i2,          i0,          i1 + 1,          i3,		s2 - 1,      i0,          i1,		i3 );//Hard constraints
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1]) >= 0 )
						m_graph.add_arc( i2,          i0,          i1 + 1,          i3,		i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1],      i0,          i1,		i3 );//Hard constraints
					else
						m_graph.add_arc( i2,          i0,          i1 + 1,          i3,		0,      i0,     i1,		i3 );//Hard constraints
				}
				//End boundary condition		
			} //for i0
		}// for i1
    }// for if
     cout<<"Finish dir-1"<<endl;
	
    //Free memory
		//vce_para.clear();
		//arc_cof.clear();
   
} //build_vce_x
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::build_vce_arcs_y(const shape_vce_type& shape_vce)
{
	cout<<"Begin vce arcs building"<<endl;
    int i0, i1, i2, i3, s0, s1, s2, s3, ii;
    int convexPower;
	//int graph_id;
	
	i3 = shape_vce.k;
	convexPower = pow_vce;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_2();
    s2 = (int)m_graph.size_1();
	
    // Construct graph arcs based on the given parameters.
    const int & bounds0 = 0;
    const int & bounds1 = 0;

	//// -- Intra-column (vertical) arcs.
    for (i2 = s2 - bounds1 - 1;
             i2 > bounds0;
             --i2) {
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
                    m_graph.add_arc(i0,          i2,          i1,          i3,          i0,          i2 - 1,          i1,          i3);
            } // for i0 
        } // for i1
    } // for i2

    // Inter-column arcs (dir-0).

	if ( i0 > 1)
	{
        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0 - 1; ++i0) {
				
				int upBound = 0, lowBound = 0; 
				//Forward
                if ( shape_vce.meanDir0[i0][i1] > 0 )
				upBound = shape_vce.meanDir0[i0][i1];
				//cout<<"The upBound is: "<<upBound<<endl;

				if ( (shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]) < 0 )
				lowBound = - ( shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]);

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.lowDir0[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.fwdCofDir0[i0][i1] * arc_weight(k),     i0,          i2,          i1,          i3,		i0 + 1,      i2 + shape_vce.meanDir0[i0][i1] - k,          i1,		i3 );
					
					m_graph.add_arc(i0,          i2,          i1,          i3,		i0 + 1,      i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1],	i1,		i3 );//Hard constraint
				
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
				{
					if ( (i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc(i0,          i2,          i1,          i3,		i0 + 1,      i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1],	i1,		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i2,          i1,          i3,		i0 + 1,      s2 - 1,	i1,		i3 );//Hard constraint
				
					
				
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1]) >= 0 )
						m_graph.add_arc(i0,          i2,          i1,          i3,		i0 + 1,      i2 + shape_vce.meanDir0[i0][i1] - shape_vce.lowDir0[i0][i1],	i1,		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i2,          i1,          i3,		i0 + 1,      0,          i1,		i3 );//Hard constraint
				}
				//End boundary condition
				//Backward
				upBound = 0, lowBound = 0;
				if ( shape_vce.meanDir0[i0][i1] < 0 )
				upBound = -shape_vce.meanDir0[i0][i1];

				if ( (shape_vce.meanDir0[i0][i1] + shape_vce.upDir0[i0][i1]) > 0 )
				lowBound = shape_vce.meanDir0[i0][i1] + shape_vce.upDir0[i0][i1];

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.upDir0[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.backCofDir0[i0][i1] * arc_weight(k),     i0 + 1,          i2,          i1,          i3,		i0,      i2 - shape_vce.meanDir0[i0][i1] - k,	i1,		i3 );

					m_graph.add_arc( i0 + 1,          i2,          i1,          i3,		i0,      i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1],	i1,		i3 );//Hard constraints
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
					{
					if ( (i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc( i0 + 1,          i2,          i1,          i3,		i0,      i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1],	i1,		i3 );//Hard constraints
					else
						m_graph.add_arc( i0 + 1,          i2,          i1,          i3,		i0,      s2 - 1,	i1,		i3 );//Hard constraints
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1]) >= 0 )
						m_graph.add_arc( i0 + 1,          i2,          i1,          i3,		i0,      i2 - shape_vce.meanDir0[i0][i1] - shape_vce.upDir0[i0][i1],	i1,		i3 );//Hard constraints
					else
						m_graph.add_arc( i0 + 1,          i2,          i1,          i3,		i0,      0,          i1,		i3 );//Hard constraints
				}
				//End boundary condition		
			} //for i0
		}// for i1
	} // for if
   cout<<"Finish dir-0"<<endl;
    // Inter-column arcs (dir-1). For 2-D case, no need
   if ( i1 > 1)
	{
    for (i1 = 0; i1 < s1 - 1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
				
				int upBound = 0, lowBound = 0; 
				//Forward
                if ( shape_vce.meanDir1[i0][i1] > 0 )
				upBound = shape_vce.meanDir1[i0][i1];

				if ( (shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]) < 0 )
				lowBound = - ( shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]);

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.lowDir1[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.fwdCofDir1[i0][i1] * arc_weight(k),     i0,          i2,          i1,          i3,		i0,     i2 + shape_vce.meanDir1[i0][i1] - k,	i1 + 1,		i3 );
					
					m_graph.add_arc( i0,          i2,          i1,          i3,		i0,      i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1], 	i1 + 1,		i3 );//Hard constraint
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
				{
					if ( (i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc(i0,          i2,          i1,          i3,		i0,          i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1],	i1 + 1,		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i2,          i1,          i3,		i0,      s2 - 1,	i1 + 1,		i3 );//Hard constraint
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1]) >= 0 )
						m_graph.add_arc(i0,          i2,          i1,          i3,		i0,       i2 + shape_vce.meanDir1[i0][i1] - shape_vce.lowDir1[i0][i1],	i1 + 1,		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i2,          i1,          i3,		i0,     0,	i1 + 1,		i3 );//Hard constraint
				}
				//End boundary condition

				//Backward
				upBound = 0, lowBound = 0;
				if ( shape_vce.meanDir1[i0][i1] < 0 )
				upBound = -shape_vce.meanDir1[i0][i1];

				if ( (shape_vce.meanDir1[i0][i1] + shape_vce.upDir1[i0][i1]) > 0 )
				lowBound = shape_vce.meanDir1[i0][i1] + shape_vce.upDir1[i0][i1];

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < shape_vce.upDir1[i0][i1]; k++)
				    m_graph.add_arc_cost(shape_vce.backCofDir1[i0][i1] * arc_weight(k),     i0,       i2,	i1 + 1,          i3,		i0,      i2 - shape_vce.meanDir1[i0][i1] - k,	i1,		i3 );
				    
					m_graph.add_arc( i0,          i2,          i1 + 1,          i3,		i0,      i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1],	i1,		i3 ); //Hard constraint
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
					{
					if ( (i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc( i0,       i2,  	i1 + 1,          i3,		i0,      i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1],	i1,		i3 );//Hard constraints
					else
						m_graph.add_arc( i0,        i2,		i1 + 1,          i3,		i0,      s2 - 1,	i1,		i3 );//Hard constraints
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1]) >= 0 )
						m_graph.add_arc( i0,	i2,		i1 + 1,		i3,		i0,      i2 - shape_vce.meanDir1[i0][i1] - shape_vce.upDir1[i0][i1],	i1,		i3 );//Hard constraints
					else
						m_graph.add_arc( i0,	i2,		i1 + 1,          i3,		i0,      0,          i1,		i3 );//Hard constraints
				}
				//End boundary condition		
			} //for i0
		}// for i1
	}//for if
    
     cout<<"Finish dir-1"<<endl;
    //Free memory
		//vce_para.clear();
		//arc_cof.clear();
   
} //build vce_y
/////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::build_graphcut_arcs()
{
    int i0, i1, i2, i3;
    int s0, s1, s2, s3;
	//int coef = 1000000;
	//int coef = 1;
	float theta = 1;   //For lymph nodes
	//float theta = 1;
	//int coef = 1;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_1();
    s2 = (int)m_graph.size_2();
    s3 = (int)m_graph.size_3();

	//Construct arcs for each node (regional term)
	for ( i3 = m_num_surf_graphsearch; i3 < s3; ++i3 ) 
		for ( i2 = 0; i2 < s2; ++i2 )
			for ( i1 = 0; i1 < s1; ++i1 )
				for ( i0 = 0; i0 < s0; ++i0 )
				{
					capacity_type cap_ob = ( capacity_type )( *m_pcost_ob )( i0, i1, i2, i3 - m_num_surf_graphsearch );
					capacity_type cap_bg = ( capacity_type )( *m_pcost_bg )( i0, i1, i2, i3 - m_num_surf_graphsearch );
					m_graph.add_st_arc( cap_ob, cap_bg, i0, i1, i2, i3);
				}

	//Construct arcs for each pair of neighboring nodes (boundary term)
    
	for ( i3 = m_num_surf_graphsearch; i3 < s3; ++i3 ) 
		for ( i2 = 1; i2 < s2 - 1; ++i2 )
			for ( i1 = 1; i1 < s1 - 1; ++i1 )
				for ( i0 = 1; i0 < s0 - 1; ++i0 )
				{
					capacity_type coef = m_neigh_coef[i3 - m_num_surf_graphsearch];
					capacity_type cap_center = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2, i3 - m_num_surf_graphsearch );
					capacity_type cap_last, cap_next;

					// (dir-0)
					cap_last = ( capacity_type )( *m_pcost_neigh )( i0 - 1, i1, i2, i3 - m_num_surf_graphsearch );
					cap_next = ( capacity_type )( *m_pcost_neigh )( i0 + 1, i1, i2, i3 - m_num_surf_graphsearch );
					m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0 - 1,      i1,          i2,		i3 );
					m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0 + 1,      i1,          i2,		i3 );
					//m_graph.add_arc_cost( coef * ( - log( 1-exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( 0.5 * 0.5 ) ) ) ),     i0,          i1,          i2,          i3,		 i0 - 1,      i1,          i2,		i3 );
					//m_graph.add_arc_cost( coef * ( - log ( 1 - exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( 0.5 * 0.5 ) ) ) ),     i0,          i1,          i2,          i3,		 i0 + 1,      i1,          i2,		i3 );

					// (dir-1)
					cap_last = ( capacity_type )( *m_pcost_neigh )( i0, i1 - 1, i2, i3 - m_num_surf_graphsearch );
					cap_next = ( capacity_type )( *m_pcost_neigh )( i0, i1 + 1, i2, i3 - m_num_surf_graphsearch );
					m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1 - 1,          i2,		i3 );
					m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1 + 1,          i2,		i3 );
					//m_graph.add_arc_cost( coef * ( - log ( 1- exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( 0.5 * 0.5 ) ) ) ),     i0,          i1,          i2,          i3,		 i0,      i1 - 1,          i2,		i3 );
					//m_graph.add_arc_cost( coef * ( - log ( 1- exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( 0.5 * 0.5 ) ) ) ),     i0,          i1,          i2,          i3,		 i0,      i1 + 1,          i2,		i3 );

					// (dir-2)
					cap_last = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2 - 1, i3 - m_num_surf_graphsearch );
					cap_next = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2 + 1, i3 - m_num_surf_graphsearch );
					m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 - 1,		i3 );
					m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 + 1 ,		i3 );
					//m_graph.add_arc_cost( coef * ( - log ( 1- exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( 0.5 * 0.5 ) ) ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 - 1,		i3 );
					//m_graph.add_arc_cost( coef * ( - log ( 1- exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( 0.5 * 0.5 ) ) ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 + 1 ,		i3 );


				}

		//Boundary condition
		//i0
	for ( i3 = m_num_surf_graphsearch; i3 < s3; ++i3 )
		for ( i2 = 1; i2 < s2 - 1; ++i2 )
			for ( i1 = 1; i1 < s1 - 1; ++i1)
			{		
					capacity_type cap_center, cap_last, cap_next;
					capacity_type coef = m_neigh_coef[i3 - m_num_surf_graphsearch];

					int i0_boundary[2];
					i0_boundary[0] = 0;
					i0_boundary[1] = s0 - 1;

					for ( int iter = 0; iter < 2; iter++ )
					{
						i0 = i0_boundary[ iter ];
						cap_center = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2, i3 - m_num_surf_graphsearch );


					    // (dir-0)
						if ( i0 - 1 >= 0 )
					    {
					       cap_last = ( capacity_type )( *m_pcost_neigh )( i0 - 1, i1, i2, i3 - m_num_surf_graphsearch );
					       m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0 - 1,     i1,          i2,		i3 );
					    }

					    if ( i0 + 1 < s0 )
					    {
					       cap_next = ( capacity_type )( *m_pcost_neigh )( i0 + 1, i1, i2, i3 - m_num_surf_graphsearch );
					       m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0 + 1,      i1,          i2,		i3 );
					    }

						// (dir-1)
					    cap_last = ( capacity_type )( *m_pcost_neigh )( i0, i1 - 1, i2, i3 - m_num_surf_graphsearch );
					    cap_next = ( capacity_type )( *m_pcost_neigh )( i0, i1 + 1, i2, i3 - m_num_surf_graphsearch );
					    m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1 - 1,          i2,		i3 );
					    m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1 + 1,          i2,		i3 );

						// (dir-2)
					    cap_last = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2 - 1, i3 - m_num_surf_graphsearch );
					    cap_next = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2 + 1, i3 - m_num_surf_graphsearch );
					    m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 - 1,		i3 );
					    m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 + 1,		i3 );

					}

			}
		//i1
	for ( i3 = m_num_surf_graphsearch; i3 < s3; ++i3 )
		for ( i2 = 1; i2 < s2 - 1; ++i2 )
			for ( i0 = 1; i0 < s0 - 1; ++i0)
			{
				
					capacity_type cap_center, cap_last, cap_next;
					capacity_type coef = m_neigh_coef[i3 - m_num_surf_graphsearch];

					int i1_boundary[2];
					i1_boundary[0] = 0;
					i1_boundary[1] = s1 - 1;

					for ( int iter = 0; iter < 2; iter++ )
					{
						i1 = i1_boundary[ iter ];
						cap_center = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2, i3 - m_num_surf_graphsearch );

						// (dir-0)
					    cap_last = ( capacity_type )( *m_pcost_neigh )( i0 - 1, i1, i2, i3 - m_num_surf_graphsearch );
					    cap_next = ( capacity_type )( *m_pcost_neigh )( i0 + 1, i1, i2, i3 - m_num_surf_graphsearch );
					    m_graph.add_arc_cost(  coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0 - 1,      i1,          i2,		i3 );
					    m_graph.add_arc_cost(  coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0 + 1,      i1,          i2,		i3 );

					    // (dir-1)
						if ( i1 - 1 >= 0 )
					    {
					       cap_last = ( capacity_type )( *m_pcost_neigh )( i0, i1 - 1, i2, i3 - m_num_surf_graphsearch );
					       m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1 - 1,          i2,		i3 );
					    }

					    if ( i1 + 1 < s1 )
					    {
					       cap_next = ( capacity_type )( *m_pcost_neigh )( i0, i1 + 1, i2, i3 - m_num_surf_graphsearch );
					       m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1 + 1,          i2,		i3 );
					    }

						// (dir-2)
					    cap_last = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2 - 1, i3 - m_num_surf_graphsearch );
					    cap_next = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2 + 1, i3 - m_num_surf_graphsearch );
					    m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 - 1,		i3 );
					    m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 + 1,		i3 );

					}

			}
		//i2
	for ( i3 = m_num_surf_graphsearch; i3 < s3; ++i3 )
		for ( i1 = 1; i1 < s1 - 1; ++i1 )
			for ( i0 = 1; i0 < s0 - 1; ++i0)
			{
				
					capacity_type cap_center, cap_last, cap_next;
					capacity_type coef = m_neigh_coef[i3 - m_num_surf_graphsearch];

					int i2_boundary[2];
					i2_boundary[0] = 0;
					i2_boundary[1] = s2 - 1;

					for ( int iter = 0; iter < 2; iter++ )
					{
						i2 = i2_boundary[ iter ];
						cap_center = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2, i3 - m_num_surf_graphsearch );
						// (dir-0)
					    cap_last = ( capacity_type )( *m_pcost_neigh )( i0 - 1, i1, i2, i3 - m_num_surf_graphsearch );
					    cap_next = ( capacity_type )( *m_pcost_neigh )( i0 + 1, i1, i2, i3 - m_num_surf_graphsearch );
					    m_graph.add_arc_cost(  coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0 - 1,      i1,          i2,		i3 );
					    m_graph.add_arc_cost(  coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0 + 1,      i1,          i2,		i3 );

					    // (dir-1)
					    cap_last = ( capacity_type )( *m_pcost_neigh )( i0, i1 - 1, i2, i3 - m_num_surf_graphsearch );
					    cap_next = ( capacity_type )( *m_pcost_neigh )( i0, i1 + 1, i2, i3 - m_num_surf_graphsearch );
					    m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1 - 1,          i2,		i3 );
					    m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1 + 1,          i2,		i3 );

					    // (dir-2)
						if ( i2 - 1 >= 0 )
					    {
					       cap_last = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2 - 1, i3 - m_num_surf_graphsearch );
					       m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_last ) * ( cap_center - cap_last ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 - 1,		i3 );
					    }

					    if ( i2 + 1 < s2 )
					    {
					       cap_next = ( capacity_type )( *m_pcost_neigh )( i0, i1, i2 + 1, i3 - m_num_surf_graphsearch );
					       m_graph.add_arc_cost( coef * exp( -1 * 0.5 * ( cap_center - cap_next ) * ( cap_center - cap_next ) / ( theta * theta ) ),     i0,          i1,          i2,          i3,		 i0,      i1,          i2 + 1,		i3 );
					    }

					}

			}

}
///////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////

/////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::build_gs_gc_arcs()
{
    int ii, i0, i1, i2, i3;
    int s0, s1, s2;

	int dir;
	for (i3 = 0; i3 < (int)m_inter_cutsearch.size(); ++i3) 
	{

        const size_type& k0 = m_inter_cutsearch[i3].k[0];
        const size_type& k1 = m_inter_cutsearch[i3].k[1];
        const int&       r = m_inter_cutsearch[i3].r;
		
		//Search for corresponding surfaces
		for (int s_index = 0; s_index < m_shape_prior.size(); s_index++  )
		{
			if (	m_shape_prior[s_index].k == k1 )
			{	
				dir = m_shape_prior[s_index].dir;
				break;
			}
		}
		
		if ( dir == 0 || dir == 1)
		{
			s0 = (int)m_graph.size_1();
			s1 = (int)m_graph.size_2();
			s2 = (int)m_graph.size_0();
			for (i1 = 0; i1 < s1; ++i1) 
				for (i0 = 0; i0 < s0; ++i0) 
					for (i2 = 0; i2 < s2; ++i2) 
					{   
						if (dir == 0)
							ii = i2 + r;
						else   //Reverse the graph
							ii = s2 - i2 + 1 + r; 
					
						if ( ii >= 0 && ii < s2 ) 
							m_graph.add_arc(i2, i0, i1, k0, ii, i0, i1, k1);
						else if ( ii >= s2 ) //Boundary case
							m_graph.add_arc(i2, i0, i1, k0, s2 - 1, i0, i1, k1);
						else //Boundary case
							m_graph.add_arc(i2, i0, i1, k0, 0, i0, i1, k1 );	
					} // for i2
		}
		else if ( dir == 2 || dir == 3)
		{
			s0 = (int)m_graph.size_0();
			s1 = (int)m_graph.size_2();
			s2 = (int)m_graph.size_1();
			for (i1 = 0; i1 < s1; ++i1) 
				for (i0 = 0; i0 < s0; ++i0) 
					for (i2 = 0; i2 < s2; ++i2) 
					{   
						if (dir == 2)
							ii = i2 + r;
						else   //Reverse the graph
							ii = s2 - i2 + 1 + r; 
					
						if ( ii >= 0 && ii < s2 ) 
							m_graph.add_arc(i0, i2, i1, k0, i0, ii, i1, k1);
						else if ( ii >= s2 ) //Boundary case
							m_graph.add_arc(i0, i2, i1, k0, i0, s2 - 1, i1, k1);
						else //Boundary case
							m_graph.add_arc(i0, i2, i1, k0, i0, 0, i1, k1 );	
					} // for i2
		}
		else if ( dir == 4 || dir == 5)
		{
			s0 = (int)m_graph.size_0();
			s1 = (int)m_graph.size_1();
			s2 = (int)m_graph.size_2();
			for (i1 = 0; i1 < s1; ++i1) 
				for (i0 = 0; i0 < s0; ++i0) 
					for (i2 = 0; i2 < s2; ++i2) 
					{   
						if (dir == 4)
							ii = i2 + r;
						else   //Reverse the graph
							ii = s2 - i2 + 1 + r; 
					
						if ( ii >= 0 && ii < s2 ) 
							m_graph.add_arc(i0, i1, i2, k0, i0, i1, ii, k1);
						else if ( ii >= s2 ) //Boundary case
							m_graph.add_arc(i0, i1, i2, k0, i0, i1, s2 - 1, k1);
						else //Boundary case
							m_graph.add_arc(i0, i1, i2, k0, i0, i1, 0, k1 );	
					} // for i2
		}
	} // for i3
} // build_gs_gc_arcs

/////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_gs_gt_multi_dir<_Cost, _Cap, _Tg>::build_gc_gc_arcs()
{
    int i0, i1, i2, i3;
    int s0, s1, s2;

	s0 = (int)m_graph.size_0();
	s1 = (int)m_graph.size_1();
	s2 = (int)m_graph.size_2();
	
	
	for (i3 = 0; i3 < (int)m_inter_cutcut.size(); ++i3) 
	{
        const size_type& k0 = m_inter_cutcut[i3].k[0];
        const size_type& k1 = m_inter_cutcut[i3].k[1];
        //const int&       r = m_inter_cutsearch[i3].r;
		for (i1 = 0; i1 < s1; ++i1) 
				for (i0 = 0; i0 < s0; ++i0) 
					for (i2 = 0; i2 < s2; ++i2) 
					{   
						capacity_type cost_0 = (*m_inter_cutcut[i3].cost_context_cut)(i0, i1, i2, 0);
						capacity_type cost_1 = (*m_inter_cutcut[i3].cost_context_cut)(i0, i1, i2, 1);
						m_graph.add_arc_cost( cost_0, i0, i1, i2, k0, i0, i1, i2, k1);				
						m_graph.add_arc_cost( cost_1, i0, i1, i2, k1, i0, i1, i2, k0);
					} // for i2
		
		
	
	} // for i3
} // build_gc_gc_arcs

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace

#endif
