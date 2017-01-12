/*
 ==========================================================================
 |   Written by Qi Song <qi-song@uiowa.edu>
 |   Department of Electrical and Computer Engineering
 |   University of Iowa
 |   
 ==========================================================================
 */

#ifndef ___OPTNET_vce_3d_terrain_MULTI_CXX___
#   define ___OPTNET_vce_3d_terrain_MULTI_CXX___

#   include <optnet/_base/except.hxx>
#   include <optnet_graphcut/optnet_vce_graphcut_terrain_multi.hxx>

#   if defined(_MSC_VER) && (_MSC_VER > 1000) && (_MSC_VER <= 1200)
#       pragma warning(disable: 4018)
#       pragma warning(disable: 4146)
#   endif
#   include <deque>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::optnet_vce_graphcut_terrain_multi() :
    m_pcost_gs(0)
{}

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::create( size_type s_0, size_type s_1, size_type s_2, size_type num_surf_graphsearch, size_type num_surf_graphcut )
{
	size_type s_3 = num_surf_graphsearch + num_surf_graphcut;

    if ( !m_graph.create( s_0, s_1, s_2, s_3 ) ) 
	{
        throw_exception(std::runtime_error(
            "optnet_vce_graphcut_terrain_multi::create: Could not create graph."
        ));
    }

    m_num_surf_graphsearch = num_surf_graphsearch;
	m_num_surf_graphcut = num_surf_graphcut;

    m_intra.resize( num_surf_graphsearch );

    // Set up default parameters.
    for ( size_type i = 0; i < num_surf_graphsearch; ++i) 
	{
        m_intra[i].tonode[0].clear();
        m_intra[i].tonode[1].clear();
        m_intra[i].dropped   = false;
        m_intra[i].climbed   = false;
        m_intra[i].circle[0] = false;
        m_intra[i].circle[1] = false;
        m_intra[i].margin[0] = 0;
        m_intra[i].margin[1] = 0;
        m_intra[i].smooth[0] = 1;
        m_intra[i].smooth[1] = 1;
    }

    m_inter.clear();

	pow_vce = 1;

	arc_cost_pos = new capacity_type* [ s_0 ];
	arc_cost_neg = new capacity_type* [ s_0 ];
	//image =new Color*[yres];
	for ( size_type i = 0; i < s_0; i++ )
	{
		arc_cost_pos[i] = new capacity_type [ s_1 ];
		arc_cost_neg[i] = new capacity_type [ s_1 ];
	}
	
    // Temporarily save a pointer to the cost vector.
    //m_pcost = &cost;

	//Initialize hard smoothness constraints;
	smooth_array.resize( num_surf_graphsearch );
	for (size_type i = 0; i < num_surf_graphsearch; ++i) {
        smooth_array[i].smoothDir0 = new int* [ s_0 ];
		smooth_array[i].smoothDir1 = new int* [ s_0 ];
		for ( size_type j = 0; j < s_0; j++)
		{
			smooth_array[i].smoothDir0[j] = new int [ s_1 ];
			smooth_array[i].smoothDir1[j] = new int [ s_1 ];
			for ( size_type k = 0; k < s_1; k++)
			{
				smooth_array[i].smoothDir0[j][k] = 1;
				smooth_array[i].smoothDir1[j][k] = 1;
			} //for k
		} //for j

    }//for i
	//Initialize vce constraints;
	vce_para.resize( num_surf_graphsearch );
	for ( size_type i = 0; i < num_surf_graphsearch; ++i ) {
        vce_para[i].meanDir0 = new int* [s_0];
		vce_para[i].meanDir1 = new int* [s_0];
		vce_para[i].upDir0 = new int* [s_0];
		vce_para[i].upDir1 = new int* [s_0];
		vce_para[i].lowDir0 = new int* [s_0];
		vce_para[i].lowDir1 = new int* [s_0];
		for ( size_type j = 0; j < s_0; j++)
		{
			vce_para[i].meanDir0[j] = new int [s_1];
		    vce_para[i].meanDir1[j] = new int [s_1];
		    vce_para[i].upDir0[j] = new int [s_1];
		    vce_para[i].upDir1[j] = new int [s_1];
		    vce_para[i].lowDir0[j] = new int [s_1];
		    vce_para[i].lowDir1[j] = new int [s_1];
			for ( size_type k = 0; k < s_1; k++)
			{
                vce_para[i].meanDir0[j][k] = 0;
		        vce_para[i].meanDir1[j][k] = 0;
		        vce_para[i].upDir0[j][k] = 1;
		        vce_para[i].upDir1[j][k] = 1;
		        vce_para[i].lowDir0[j][k] = 1;
		        vce_para[i].lowDir1[j][k] = 1;
			} //for k
		} //for j

    }//for i
	//Initialize arc coefficients;
	arc_cof.resize( num_surf_graphsearch );
	for ( size_type i = 0; i < num_surf_graphsearch; ++i) {
        arc_cof[i].fwdDir0 = new int* [s_0];
		arc_cof[i].fwdDir1 = new int* [s_0];
		arc_cof[i].backDir0 = new int* [s_0];
		arc_cof[i].backDir1 = new int* [s_0];

		for ( size_type j = 0; j < s_0; j++)
		{
			arc_cof[i].fwdDir0[j] = new int [s_1];
		    arc_cof[i].fwdDir1[j] = new int [s_1];
		    arc_cof[i].backDir0[j] = new int [s_1];
		    arc_cof[i].backDir1[j] = new int [s_1];

			for ( size_type k = 0; k < s_1; k++)
			{
                arc_cof[i].fwdDir0[j][k] = 0;
		        arc_cof[i].fwdDir1[j][k] = 0;
		        arc_cof[i].backDir0[j][k] = 0;
		        arc_cof[i].backDir1[j][k] = 0;

			} //for k
		} //for j

    }//for i
}
///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
int
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::arc_weight(int k , size_type i0, size_type i1, size_type i3, int dir, int fwdFlag )
                                          
{
   int weight;
   int cof;
   //power = pow_vce;
   if ( dir == 0 && fwdFlag == 0)
   cof = arc_cof[i3].fwdDir0[i0][i1];
   else if ( dir == 0 && fwdFlag == 1)
   cof = arc_cof[i3].backDir0[i0][i1];
   else if ( dir == 1 && fwdFlag == 0)
   cof = arc_cof[i3].fwdDir1[i0][i1];
   else
   cof = arc_cof[i3].backDir1[i0][i1];

   if ( k == 0)
   weight = cof;
   else
   weight = cof * ( pow( (k+1), pow_vce )-2 * pow( k, pow_vce) + pow( (k-1), pow_vce) );
   return weight;
}


//////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::set_cutsearch_relation(size_type k0,
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

///////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::solve_all(net_base_type& net, 
                                            capacity_type* pflow
                                            )
{
    size_type       i0, i1, i2, i3;
    capacity_type   flow;

    if (net.size_0() != m_graph.size_0() || 
        net.size_1() != m_graph.size_1() ||
		net.size_2() != m_graph.size_2() ||
        net.size_3() != m_graph.size_3()
        ) {
        // Throw an invalid_argument exception.
        throw_exception(
            std::invalid_argument(
            "optnet_vce_graphcut_terrain_multi::solve: The output image size must match the graph size."
        ));
    }

	// Assign the cost of graph nodes based on the input cost
    // vector. We also perform the "translation operation"
    // here to guaranttee a non-empty solution.
    transform_costs();
    cout<<"Finish cost transform"<<endl;

    // Build the arcs of the graphs.
	build_vce_arcs();
    build_graphcut_arcs();
	build_intersurface_arcs();
	cout << "Finish build arcs" << endl;

    // Calculate max-flow/min-cut.
    flow = m_graph.solve();
    
    // Get the labeled image..
    for (i3 = 0; i3 < m_graph.size_3(); ++i3) {	
            for (i1 = 0; i1 < m_graph.size_1(); ++i1) {
                for (i0 = 0; i0 < m_graph.size_0(); ++i0) {

                    // Find upper envelope.
					if ( i3 < m_num_surf_graphsearch )
					{
						 size_type lowest;
                
                         lowest = m_intra[i3].margin[0];

                         for ( i2 = lowest; i2 < m_graph.size_2() - m_intra[i3].margin[1]; ++i2) 
						 {
                    
                             // Find upper envelope.
                             if (!m_graph.in_source_set(i0, i1, i2, i3))
                             break;
                         }
						 if ( i3 == 0)
							 net( i0, i1, i2, i3 ) = 1;
						 else 
							 net( i0, i1, m_graph.size_2() - i2 - 1, i3 ) = 1;
					}
					else
					{
                      for (i2 = 0; i2 < m_graph.size_2(); ++i2) 
					  {
						
                         if ( m_graph.in_source_set( i0, i1, i2, i3 ) )
						    net( i0, i1, i2, i3 ) = 1;
					     else
						    net( i0, i1, i2, i3 ) = 0;
							
							
					  }
					}

                } // for i0
            } // for i1
    } // for i3

    if (0 != pflow)
        *pflow = flow;
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::get_bounds_of_subgraphs()
{
    bool                stop[2];
    int                 r, new_bound;
    size_type           i, k, n = m_intra.size();
    std::deque<_Intra*> Q[2];

    //
    // Hopefully this is also a correct way of detecting illegal surface
    // relation specifications. 
    //

    //
    // STEP 1: Initialize.
    //
    for (i = 0; i < n; ++i) {

        _Intra* p = &m_intra[i];
        
        if (p->climbed) p->visits[0] = 0;
        else { // else 1
            Q[0].push_back(p);
            p->visits[0] = 1;
        } // else 1

        if (p->dropped) p->visits[1] = 0;
        else { // else 2
            Q[1].push_back(p);
            p->visits[1] = 1;
        } // else 2

    } // for

    //
    // STEP 2: Compute the margin of the graphs.
    //
    if (Q[0].empty() || Q[1].empty()) {
        // The specified relations are invalid.
        throw_exception(
            std::invalid_argument(
            "optnet_vce_graphcut_terrain_multi::get_bounds_of_subgraphs: Invalid surface interrelation."
        ));
    }
    else {
        stop[0] = false;
        stop[1] = false;
        
        //
        // Lower margin.
        do {
            _Intra* p = Q[0].front();
            Q[0].pop_front();

            for (i = 0; i < p->tonode[0].size(); ++i) {
                k = p->tonode[0][i].first;
                r = p->tonode[0][i].second;
                new_bound = p->margin[0] + r;
                if (m_intra[k].margin[0] < new_bound) {
                    if (++m_intra[k].visits[0] < (int)n) {
                        m_intra[k].margin[0] = new_bound;
                        Q[0].push_back(&m_intra[k]);
                    }
                    else {
                        stop[0] = true;
                        break;
                    } // else
                } // if
            } // for

        } // do
        while (!Q[0].empty() && !stop[0]);

        //
        // Upper margin.
        do {
            _Intra* p = Q[1].front();
            Q[1].pop_front();

            for (i = 0; i < p->tonode[1].size(); ++i) {
                k = p->tonode[1][i].first;
                r = p->tonode[1][i].second;
                new_bound = p->margin[1] + r;
                if (m_intra[k].margin[1] < new_bound) {
                    if (++m_intra[k].visits[1] < (int)n) {
                        m_intra[k].margin[1] = new_bound;
                        Q[1].push_back(&m_intra[k]);
                    }
                    else {
                        stop[1] = true;
                        break;
                    } // else
                } // if
            } // for

        } // do
        while (!Q[1].empty() && !stop[1]);

        if (stop[0] || stop[1]) {
            // There is a loop in the relations specified.
            throw_exception(
                std::invalid_argument(
                "optnet_vce_graphcut_terrain_multi::get_bounds_of_subgraphs: Invalid surface interrelation."
            ));
        } // if
    }
}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::transform_costs()
{
    assert(m_pcost_gs != 0);

    size_type            i0, i1, i2, i3;

    const size_type& s0 = m_pcost_gs->size_0();
    const size_type& s1 = m_pcost_gs->size_1();
    const size_type& s2 = m_pcost_gs->size_2();
    const size_type& s3 = m_pcost_gs->size_3();

    m_graph.set_initial_flow(0);

        // Construct the s-t graph "G_st".
        for (i3 = 0; i3 < s3; ++i3) {
            // Nodes below margin[0] and above s2-margin[1] will
            // be ignored.
            for (i2 = m_intra[i3].margin[0] + 1;
                i2 < s2 - m_intra[i3].margin[1];
                ++i2) {
                for (i1 = 0; i1 < s1; ++i1) {
                    for (i0 = 0; i0 < s0; ++i0) {
						capacity_type cap;
                   
						if ( i3 == 0 )
                        cap = (capacity_type)(*m_pcost_gs)(i0, i1, i2, i3)
                            - (capacity_type)(*m_pcost_gs)(i0, i1, i2 - 1, i3);
						else     //Reverse the image. 
						cap = (capacity_type)(*m_pcost_gs)(i0, i1, s2 - i2 - 1, i3)
                            - (capacity_type)(*m_pcost_gs)(i0, i1, s2 - i2, i3);

                        if (cap >= 0) // non-negative -> connect to t
                            m_graph.add_st_arc(0, +cap, i0, i1, i2, i3);
                        else          // negative     -> connect to s
                            m_graph.add_st_arc(-cap, 0, i0, i1, i2, i3);

                    } // for i0
                } // for i1
            } // for i2
        } // for i3

        for (i3 = 0; i3 < s3; ++i3) {
            i2 = m_intra[i3].margin[0]; // lowest margin
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {
                    // cost = -1
                    m_graph.add_st_arc(100000, 0, i0, i1, i2, i3);
                } // for i0
            } // for i1
        } // for i3



}

///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::build_arcs()
{
    int ii, i0, i1, i2, i3, i2a, i2b;
    int s0, s1, s2, s3;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_1();
    s2 = (int)m_graph.size_2();
    s3 = (int)(m_num_surf_graphsearch);

    // Clear all constructed arcs.
    //m_graph.clear_arcs();

    //
    // Construct intra-surface arcs here.
    //

    for (i3 = 0; i3 < s3; ++i3) {

        const bool& circle0 = m_intra[i3].circle[0];
        const bool& circle1 = m_intra[i3].circle[1];
        const int & bounds0 = m_intra[i3].margin[0];
        const int & bounds1 = m_intra[i3].margin[1];
        const int & smooth0 = m_intra[i3].smooth[0];
        const int & smooth1 = m_intra[i3].smooth[1];

        // -- Intra-column (vertical) arcs.
        for (i2 = s2 - bounds1 - 1;
             i2 > bounds0;
             --i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 0; i0 < s0; ++i0) {
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1,          i2 - 1,          i3);
                } // for i0 
            } // for i1
        } // for i2

        // -- Inter-column arcs (dir-0).
	   if ( i0 > 1)
	   {
        
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 1; i0 < s0 - 1; ++i0) {
					for (i2 = s2 - bounds1 - 1; i2 > smooth_array[i3].smoothDir0[i0][i1] + bounds0; --i2) {
					
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 - 1,      i1,          i2 - smooth_array[i3].smoothDir0[i0][i1],    i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 + 1,      i1,          i2 - smooth_array[i3].smoothDir0[i0][i1],    i3);
					} }{
                    for (i2 = s2 - bounds1 - 1; i2 > smooth_array[i3].smoothDir0[i0][i1] + bounds0; --i2) {
                    m_graph.add_arc(0,           i1,          i2,          i3,          1,           i1,          i2 - smooth_array[i3].smoothDir0[i0][i1],    i3);
                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          s0 - 2,      i1,          i2 - smooth_array[i3].smoothDir0[i0][i1],    i3);
					}
                } // for i0 boundary condition
            } // for i1
	   }// for if 
        
	   // -- Inter-column arcs (dir-1).
       if ( i1 > 1 )
	   {
			 
            for (i0 = 0; i0 < s0; ++i0) {
                for (i1 = 1; i1 < s1 - 1; ++i1) {
					for (i2 = s2 - bounds1 - 1; i2 > smooth_array[i3].smoothDir1[i0][i1] + bounds0; --i2) {
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 - 1,      i2 - smooth_array[i3].smoothDir1[i0][i1],    i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 + 1,      i2 - smooth_array[i3].smoothDir1[i0][i1],    i3);
			        }
                   } {
                   for (i2 = s2 - bounds1 - 1;i2 > smooth_array[i3].smoothDir1[i0][i1] + bounds0;--i2) {
                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          1,           i2 - smooth_array[i3].smoothDir1[i0][i1],    i3);
                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          s1 - 2,      i2 - smooth_array[i3].smoothDir1[i0][i1],    i3);
			      }
                } // for i0
            } // for i1
	   }// for if
        // -- The base-set
			 
        { 
            for (i1 = 0; i1 < s1; ++i1) {
                for (i0 = 1; i0 < s0 - 1; ++i0) {


                    i2  = i2a = i2b = bounds0;

                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 - 1,      i1,          i2a,             i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0 + 1,      i1,          i2b,             i3);
                } {
                    // i0 = 0

                    i2  = i2b = bounds0;

                    m_graph.add_arc(0,           i1,          i2,          i3,          1,           i1,          i2b,             i3);

                    // i0 = s0 - 1

                        i2  = i2a = bounds0;

                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          s0 - 2,      i1,          i2a,             i3);
                } // for i0
            } // for i1
        }
		/*

        {
            for (i0 = 0; i0 < s0; ++i0) {
                for (i1 = 1; i1 < s1 - 1; ++i1) {
                    

                        i2  = i2a = i2b = bounds0;

                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 - 1,      i2a,             i3);
                    m_graph.add_arc(i0,          i1,          i2,          i3,          i0,          i1 + 1,      i2b,             i3);
                } {
                    // i1 = 0

                        i2  = i2b = bounds0;

                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          1,           i2b,             i3);

                    // i1 = s1 - 1

                        i2  = i2a = bounds0;

                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          s1 - 2,      i2a,             i3);
                } // for i0
            } // for i1
        }
		*/

        // -- Circular graph connections.
        if (circle0) {
            for (i2 = s2 - bounds1 - 1;
                 i2 > smooth0 + bounds0;
                 --i2) {
                for (i1 = 0; i1 < s1; ++i1) {
                    m_graph.add_arc(0,           i1,          i2,          i3,          s0 - 1,      i1,          i2 - smooth0,    i3);
                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          0,           i1,          i2 - smooth0,    i3);
                }
            } {
                for (i1 = 0; i1 < s1; ++i1) {
                    // i0 = 0

                        i2  = i2a = bounds0;

                    m_graph.add_arc(0,           i1,          i2,          i3,          s0 - 1,      i1,          i2a,             i3);
 
                    // i0 = s0 - 1

                        i2  = i2b = bounds0;

                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          0,           i1,          i2b,             i3);
                }
            }
        } // if

        if (circle1) {
            for (i2 = s2 - bounds1 - 1;
                 i2 > smooth1 + bounds0;
                 --i2) {
                for (i0 = 0; i0 < s0; ++i0) {
                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          s1 - 1,      i2 - smooth1,    i3);
                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          0,           i2 - smooth1,    i3);
                }                                                          
            } {                                      
                for (i0 = 0; i0 < s0; ++i0) {                              
                    // i1 = 0

                        i2  = i2a = bounds0;

                    m_graph.add_arc(i0,          0,           i2,          i3,          i0,          s1 - 1,      i2a,             i3);

                    // i1 = s1 - 1

                        i2  = i2b = bounds0;

                    m_graph.add_arc(i0,          s1 - 1,      i2,          i3,          i0,          0,           i2b,            i3);
                }
            }
        } // if

    } // for i3
    
    // Construct inter-surface arcs here.
	
    for (i3 = 0; i3 < (int)m_inter.size(); ++i3) {

        const size_type& k0 = m_inter[i3].k[0];
        const size_type& k1 = m_inter[i3].k[1];
        const int&       r0 = m_inter[i3].r[0];
        const int&       r1 = m_inter[i3].r[1];

        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {

                // k0 --> k1
                for (i2 = m_intra[k0].margin[0];
                     i2 < s2 - m_intra[k0].margin[1];
                     ++i2) {
                    
                    ii = i2 + r0;
                    if (ii <= m_intra[k1].margin[0] || 
                        ii >= s2 - m_intra[k1].margin[1]) continue;
                    m_graph.add_arc(i0, i1, i2, k0, i0, i1, ii, k1);
                }

                // k1 --> k0
                for (i2 = m_intra[k1].margin[0];
                     i2 < s2 - m_intra[k1].margin[1];
                     ++i2) {
                    
                    ii = i2 + r1;
                    if (ii <= m_intra[k0].margin[0] || 
                        ii >= s2 - m_intra[k0].margin[1]) continue;
                    m_graph.add_arc(i0, i1, i2, k1, i0, i1, ii, k0);
                }

            } // for i0
        } // for i1
    } // for i3
	

    // Interconnect "zero-plane"s.
	
    for (i3 = 1; i3 < s3; ++i3) {

        m_graph.add_arc(0, 0, m_intra[i3    ].margin[0], i3,     
                        0, 0, m_intra[i3 - 1].margin[0], i3 - 1
                        );
        m_graph.add_arc(0, 0, m_intra[i3 - 1].margin[0], i3 - 1,
                        0, 0, m_intra[i3    ].margin[0], i3
                        );
    }
	

}
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::build_vce_arcs()
{
	cout<<"Begin vce arcs building"<<endl;
    int i0, i1, i2, i3, s0, s1, s2, s3, ii;
    int convexPower;

	convexPower = pow_vce;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_1();
    s2 = (int)m_graph.size_2();
	s3 = (int)(m_num_surf_graphsearch);

    // Clear all constructed arcs.
    //m_graph.clear_arcs();

    //
    // Construct graph arcs based on the given parameters.
    // -- Sorry, the lines here are a bit too long. ;-)
    //
   for (i3 = 0; i3 < s3; ++i3) {

        const bool& circle0 = m_intra[i3].circle[0];
        const bool& circle1 = m_intra[i3].circle[1];
        const int & bounds0 = m_intra[i3].margin[0];
        const int & bounds1 = m_intra[i3].margin[1];
        //const int & smooth0 = m_intra[i3].smooth[0];
        //const int & smooth1 = m_intra[i3].smooth[1];
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
                if ( vce_para[i3].meanDir0[i0][i1] > 0 )
				upBound = vce_para[i3].meanDir0[i0][i1];
				//cout<<"The upBound is: "<<upBound<<endl;

				if ( (vce_para[i3].meanDir0[i0][i1] - vce_para[i3].lowDir0[i0][i1]) < 0 )
				lowBound = - ( vce_para[i3].meanDir0[i0][i1] - vce_para[i3].lowDir0[i0][i1]);

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < vce_para[i3].lowDir0[i0][i1]; k++)
				    m_graph.add_arc_cost(arc_weight(k, i0, i1, i3,  0,   0),     i0,          i1,          i2,          i3,		i0 + 1,      i1,          i2 + vce_para[i3].meanDir0[i0][i1] - k,		i3 );
					
					m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          i2 + vce_para[i3].meanDir0[i0][i1] - vce_para[i3].lowDir0[i0][i1],		i3 );//Hard constraint
				
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
				{
					if ( (i2 + vce_para[i3].meanDir0[i0][i1] - vce_para[i3].lowDir0[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          i2 + vce_para[i3].meanDir0[i0][i1] - vce_para[i3].lowDir0[i0][i1],		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          s2 - 1,		i3 );//Hard constraint
				
					
				
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 + vce_para[i3].meanDir0[i0][i1] - vce_para[i3].lowDir0[i0][i1]) >= 0 )
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          i2 + vce_para[i3].meanDir0[i0][i1] - vce_para[i3].lowDir0[i0][i1],		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0 + 1,      i1,          0,		i3 );//Hard constraint
				}
				//End boundary condition
				//Backward
				upBound = 0, lowBound = 0;
				if ( vce_para[i3].meanDir0[i0][i1] < 0 )
				upBound = -vce_para[i3].meanDir0[i0][i1];

				if ( (vce_para[i3].meanDir0[i0][i1] + vce_para[i3].upDir0[i0][i1]) > 0 )
				lowBound = vce_para[i3].meanDir0[i0][i1] + vce_para[i3].upDir0[i0][i1];

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < vce_para[i3].upDir0[i0][i1]; k++)
				    m_graph.add_arc_cost(arc_weight(k, i0, i1, i3, 0,  1),     i0 + 1,          i1,          i2,          i3,		i0,      i1,          i2 - vce_para[i3].meanDir0[i0][i1] - k,		i3 );

					m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          i2 - vce_para[i3].meanDir0[i0][i1] - vce_para[i3].upDir0[i0][i1],		i3 );//Hard constraints
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
					{
					if ( (i2 - vce_para[i3].meanDir0[i0][i1] - vce_para[i3].upDir0[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          i2 - vce_para[i3].meanDir0[i0][i1] - vce_para[i3].upDir0[i0][i1],		i3 );//Hard constraints
					else
						m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          s2 - 1,		i3 );//Hard constraints
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 - vce_para[i3].meanDir0[i0][i1] - vce_para[i3].upDir0[i0][i1]) >= 0 )
						m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          i2 - vce_para[i3].meanDir0[i0][i1] - vce_para[i3].upDir0[i0][i1],		i3 );//Hard constraints
					else
						m_graph.add_arc( i0 + 1,          i1,          i2,          i3,		i0,      i1,          0,		i3 );//Hard constraints
				}
				//End boundary condition

				
				
			} //for i0
			
		}// for i1

		//tmp for circle case
		if (circle0) {
            for (i2 = 0; i2 < s2; i2++ )
                for (i1 = 0; i1 < s1; ++i1) {
                    m_graph.add_arc(0,           i1,          i2,          i3,          s0 - 1,      i1,          i2,    i3);
                    m_graph.add_arc(s0 - 1,      i1,          i2,          i3,          0,           i1,          i2,    i3);
                }
            
		}
		///tmp
  
	  //Smaller than hard smoothness constraints
	  /*
		for (i2 = smooth_array[i3].smoothDir0[i0][i1]; i2>0; --i2){
			for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 1; i0 < s0 - 1; ++i0) {

				m_graph.add_arc_cost(arc_cost_pos[i0][i1],     i0,          i1,          i2,		i3,          i0 + 1,      i1,          i2,		i3 );
				m_graph.add_arc_cost(arc_cost_neg[i0][i1],     i0,          i1,          i2,		i3,          i0 - 1,      i1,          i2,		i3 );
				
				for (int k = 1; k < i2; k++){
				m_graph.add_arc_cost(arc_cost_pos[i0][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ), i0,          i1,          i2,          i3,		i0 + 1,      i1,          i2 - k,		i3);
				m_graph.add_arc_cost(arc_cost_neg[i0][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ), i0,          i1,          i2,          i3,		i0 - 1,      i1,          i2 - k,		i3);
                //m_graph.add_arc(i0,          i1,          i2,          i0 - 1,      i1,          i2 - m_smooth[0]);
                //m_graph.add_arc(i0,          i1,          i2,          i0 + 1,      i1,          i2 - m_smooth[0]);
				}// for k
            
			}// for i0 
			{ //boundary case
                //m_graph.add_arc(0,           i1,          i2,          1,           i1,          i2 - m_smooth[0]);
                //m_graph.add_arc(s0 - 1,      i1,          i2,          s0 - 2,      i1,          i2 - m_smooth[0]);

				m_graph.add_arc_cost(arc_cost_pos[0][i1],     0,          i1,          i2,		i3,           1,      i1,          i2,		i3 );
				m_graph.add_arc_cost(arc_cost_neg[s0-1][i1],     s0-1,          i1,          i2,		i3,        s0-2,      i1,          i2,		i3 );
				
				for (int k = 1; k < i2; k++){
				m_graph.add_arc_cost(arc_cost_pos[0][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ),     0,          i1,          i2,		i3,           1,      i1,          i2 - k,		i3);
				m_graph.add_arc_cost(arc_cost_neg[s0-1][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ),     s0-1,          i1,          i2,		i3,        s0-2,      i1,          i2 - k,		i3);
                } // for k
              }// for boundary
			}// for i1

		}// for i2
		*/
	} // for if
     
   cout<<"Finish dir-0"<<endl;
    // Inter-column arcs (dir-1). For 2-D case, no need
   if ( i1 > 1)
	{
    for (i1 = 0; i1 < s1 - 1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {
				
				int upBound = 0, lowBound = 0; 
				//Forward
                if ( vce_para[i3].meanDir1[i0][i1] > 0 )
				upBound = vce_para[i3].meanDir1[i0][i1];

				if ( (vce_para[i3].meanDir1[i0][i1] - vce_para[i3].lowDir1[i0][i1]) < 0 )
				lowBound = - ( vce_para[i3].meanDir1[i0][i1] - vce_para[i3].lowDir1[i0][i1]);

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < vce_para[i3].lowDir1[i0][i1]; k++)
				    m_graph.add_arc_cost(arc_weight(k, i0, i1, i3, 1,  0),     i0,          i1,          i2,          i3,		i0,      i1 + 1,          i2 + vce_para[i3].meanDir1[i0][i1] - k,		i3 );
					
					m_graph.add_arc( i0,          i1,          i2,          i3,		i0,      i1 + 1,          i2 + vce_para[i3].meanDir1[i0][i1] - vce_para[i3].lowDir1[i0][i1],		i3 );//Hard constraint
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
				{
					if ( (i2 + vce_para[i3].meanDir1[i0][i1] - vce_para[i3].lowDir1[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0,      i1 + 1,          i2 + vce_para[i3].meanDir1[i0][i1] - vce_para[i3].lowDir1[i0][i1],		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0,      i1 + 1,          s2 - 1,		i3 );//Hard constraint
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 + vce_para[i3].meanDir1[i0][i1] - vce_para[i3].lowDir1[i0][i1]) >= 0 )
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0,      i1 + 1,          i2 + vce_para[i3].meanDir1[i0][i1] - vce_para[i3].lowDir1[i0][i1],		i3 );//Hard constraint
					else
						m_graph.add_arc(i0,          i1,          i2,          i3,		i0,      i1 + 1,          0,		i3 );//Hard constraint
				}
				//End boundary condition

				//Backward
				upBound = 0, lowBound = 0;
				if ( vce_para[i3].meanDir1[i0][i1] < 0 )
				upBound = -vce_para[i3].meanDir1[i0][i1];

				if ( (vce_para[i3].meanDir1[i0][i1] + vce_para[i3].upDir1[i0][i1]) > 0 )
				lowBound = vce_para[i3].meanDir1[i0][i1] + vce_para[i3].upDir1[i0][i1];

				for (i2 = s2 - 1 - upBound; i2 >= lowBound; --i2) {
                
					for (int k = 0; k < vce_para[i3].upDir1[i0][i1]; k++)
				    m_graph.add_arc_cost(arc_weight(k, i0, i1, i3, 1,   1),     i0,          i1 + 1,          i2,          i3,		i0,      i1,          i2 - vce_para[i3].meanDir1[i0][i1] - k,		i3 );
				    
					m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          i2 - vce_para[i3].meanDir1[i0][i1] - vce_para[i3].upDir1[i0][i1],		i3 ); //Hard constraint
				} //for i2
				//Boundary condition
				for ( i2 = s2 - upBound; i2 <= s2 - 1; i2++ )
					{
					if ( (i2 - vce_para[i3].meanDir1[i0][i1] - vce_para[i3].upDir1[i0][i1]) <= (s2 - 1) )
						m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          i2 - vce_para[i3].meanDir1[i0][i1] - vce_para[i3].upDir1[i0][i1],		i3 );//Hard constraints
					else
						m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          s2 - 1,		i3 );//Hard constraints
				}
				for ( i2 = lowBound - 1; i2 >= 0; i2-- )
				{
					if ( (i2 - vce_para[i3].meanDir1[i0][i1] - vce_para[i3].upDir1[i0][i1]) >= 0 )
						m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          i2 - vce_para[i3].meanDir1[i0][i1] - vce_para[i3].upDir1[i0][i1],		i3 );//Hard constraints
					else
						m_graph.add_arc( i0,          i1 + 1,          i2,          i3,		i0,      i1,          0,		i3 );//Hard constraints
				}
				//End boundary condition

				
				
			} //for i0
			
		}// for i1

	//tmp for circle case
		if (circle1) {
            for (i2 = 0; i2 < s2; i2++ )
                for (i0 = 0; i0 < s0; ++i0) {
                    m_graph.add_arc(i0,           0,          i2,          i3,          i0,      s1 - 1,          i2,    i3);
                    m_graph.add_arc(i0,      s1 - 1,          i2,          i3,          i0,           0,          i2,    i3);
                }
            
		}
   
	  //Smaller than hard smoothness constraints
      /*
		for (i2 = smooth_array[i3].smoothDir1[i0][i1]; i2>0; --i2){
		  for (i0 = 0; i0 < s0; ++i0) {
			for (i1 = 1; i1 < s1 - 1; ++i1) {
            

				m_graph.add_arc_cost(arc_cost_pos[i0][i1],     i0,          i1,          i2,		i3,          i0,      i1 + 1,          i2,		i3 );
				m_graph.add_arc_cost(arc_cost_neg[i0][i1],     i0,          i1,          i2,		i3,          i0,      i1 - 1,          i2,		i3 );
				
				for (int k = 1; k < i2; k++){
				m_graph.add_arc_cost(arc_cost_pos[i0][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ), i0,          i1,          i2,          i3,		i0,      i1 + 1,          i2 - k,		i3);
				m_graph.add_arc_cost(arc_cost_neg[i0][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ), i0,          i1,          i2,          i3,		i0,      i1 - 1,          i2 - k,		i3);
                //m_graph.add_arc(i0,          i1,          i2,          i0 - 1,      i1,          i2 - m_smooth[0]);
                //m_graph.add_arc(i0,          i1,          i2,          i0 + 1,      i1,          i2 - m_smooth[0]);
				}// for k
            
			}// for i1 
			{ //boundary case
                //m_graph.add_arc(0,           i1,          i2,          1,           i1,          i2 - m_smooth[0]);
                //m_graph.add_arc(s0 - 1,      i1,          i2,          s0 - 2,      i1,          i2 - m_smooth[0]);

				m_graph.add_arc_cost(arc_cost_pos[i0][0],     i0,          0,          i2,		i3,           i0,      1,          i2,		i3 );
				m_graph.add_arc_cost(arc_cost_neg[i0][s1-1],     i0,          s1 - 1,          i2,		i3,        i0,      s1 - 2,          i2,		i3 );
				
				for (int k = 1; k < i2; k++){
				m_graph.add_arc_cost(arc_cost_pos[i0][0] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ),     i0,          0,          i2,		i3,           i0,      1,          i2 - k,		i3);
				m_graph.add_arc_cost(arc_cost_neg[i0][s1-1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ),     i0,          s1 - 1,          i2,		i3,        i0,      s1 - 2,          i2 - k,		i3);
                } // for k
              }// for boundary
			
			}// for i0
			

		}// for i2
		*/
	} // for if
    
     cout<<"Finish dir-1"<<endl;
    // Circular graph connections.
    /*
    if (circle0) {
        for (i2 = s2 - 1; i2 > smooth0; --i2) {
            for (i1 = 0; i1 < s1; ++i1) {
                //m_graph.add_arc(0,           i1,          i2,          s0 - 1,      i1,          i2 - m_smooth[0]);
                //m_graph.add_arc(s0 - 1,      i1,          i2,          0,           i1,          i2 - m_smooth[0]);
				m_graph.add_arc_cost(arc_cost_pos[s0-1][i1],     s0 - 1,          i1,          i2,		i3,          0,      i1,          i2,		i3 );
				m_graph.add_arc_cost(arc_cost_neg[0][i1],     0,          i1,          i2,         i3,		 s0 - 1,      i1,          i2,		i3 );
				
				for (int k = 1; k < smooth0; k++){
				m_graph.add_arc_cost(arc_cost_pos[s0-1][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ), s0-1,          i1,          i2,		i3,          0,      i1,          i2 - k,		i3);
				m_graph.add_arc_cost(arc_cost_neg[0][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ), 0,          i1,          i2,		i3,          s0 - 1,      i1,          i2 - k,		i3);
                //m_graph.add_arc(i0,          i1,          i2,          i0 - 1,      i1,          i2 - m_smooth[0]);
                //m_graph.add_arc(i0,          i1,          i2,          i0 + 1,      i1,          i2 - m_smooth[0]);
				}// for k
            }//for i1
        } //for i2
		
		{ 
			//Smaller than hard smoothness constraints
		    for (i2 = smooth0; i2>0; --i2){
				for (i1 =0; i1 < s1; ++i1)
				{
				m_graph.add_arc_cost(arc_cost_pos[s0-1][i1],     s0 - 1,          i1,          i2,		i3,          0,      i1,          i2,		i3 );
				m_graph.add_arc_cost(arc_cost_neg[0][i1],     0,          i1,          i2,		i3,          s0 - 1,      i1,          i2,		i3 );
		 	
				for (int k = 1; k < i2; k++){
				m_graph.add_arc_cost(arc_cost_pos[s0-1][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ), s0-1,          i1,          i2,		i3,          0,      i1,          i2 - k,		i3);
				m_graph.add_arc_cost(arc_cost_neg[0][i1] * ( pow( (k+1), convexPower )-2 * pow( k,convexPower) + pow( (k-1), convexPower) ), 0,          i1,          i2,		i3,          s0 - 1,      i1,          i2 - k,		i3);
                }//for k
			    } // for i1
			}// for i2
		}
		
    } // if circle
	*/

  } // for i3

  // Construct inter-surface arcs here.
	
    for (i3 = 0; i3 < (int)m_inter.size(); ++i3) {

        const size_type& k0 = m_inter[i3].k[0];
        const size_type& k1 = m_inter[i3].k[1];
        const int&       r0 = m_inter[i3].r[0];
        const int&       r1 = m_inter[i3].r[1];

        for (i1 = 0; i1 < s1; ++i1) {
            for (i0 = 0; i0 < s0; ++i0) {

                // k0 --> k1
                for (i2 = m_intra[k0].margin[0];
                     i2 < s2 - m_intra[k0].margin[1];
                     ++i2) {
                    
                    ii = i2 + r0;
                    if (ii <= m_intra[k1].margin[0] || 
                        ii >= s2 - m_intra[k1].margin[1]) continue;
                    m_graph.add_arc(i0, i1, i2, k0, i0, i1, ii, k1);
                }

                // k1 --> k0
                for (i2 = m_intra[k1].margin[0];
                     i2 < s2 - m_intra[k1].margin[1];
                     ++i2) {
                    
                    ii = i2 + r1;
                    if (ii <= m_intra[k0].margin[0] || 
                        ii >= s2 - m_intra[k0].margin[1]) continue;
                    m_graph.add_arc(i0, i1, i2, k1, i0, i1, ii, k0);
                }

            } // for i0
        } // for i1
    } // for i3
    
	// Interconnect "zero-plane"s.
    /*
    for (i3 = 1; i3 < s3; ++i3) {

        m_graph.add_arc(0, 0, m_intra[i3    ].margin[0], i3,     
                        0, 0, m_intra[i3 - 1].margin[0], i3 - 1
                        );
        m_graph.add_arc(0, 0, m_intra[i3 - 1].margin[0], i3 - 1,
                        0, 0, m_intra[i3    ].margin[0], i3
                        );
	  }

     */	
    //Free memory
		vce_para.clear();
		arc_cof.clear();


  
   
} //build vce

/////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////
template <typename _Cost, typename _Cap, typename _Tg>
void
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::build_graphcut_arcs()
{
    int i0, i1, i2, i3;
    int s0, s1, s2, s3;
	int coef = 1000000;
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
optnet_vce_graphcut_terrain_multi<_Cost, _Cap, _Tg>::build_intersurface_arcs()
{
    int ii, i0, i1, i2, i3;
    int s0, s1, s2;

    s0 = (int)m_graph.size_0();
    s1 = (int)m_graph.size_1();
    s2 = (int)m_graph.size_2();
    
	for (i3 = 0; i3 < (int)m_inter_cutsearch.size(); ++i3) 
	{

        const size_type& k0 = m_inter_cutsearch[i3].k[0];
        const size_type& k1 = m_inter_cutsearch[i3].k[1];
        const int&       r = m_inter_cutsearch[i3].r;
		
        for (i1 = 0; i1 < s1; ++i1) 
		{
            for (i0 = 0; i0 < s0; ++i0) 
			{
               // k0 --> k1
                for (i2 = 0; i2 < s2; ++i2) 
				{   
					if (k1 == 0)
                      ii = i2 + r;
					else if (k1 == 1)  //Reverse the graph
					  ii = s2 - i2 + 1 + r; 
					
                      if ( ii >= 0 && ii < s2 ) 
                         m_graph.add_arc(i0, i1, i2, k0, i0, i1, ii, k1);
					  else if ( ii >= s2 ) //Boundary case
					     m_graph.add_arc(i0, i1, i2, k0, i0, i1, s2 - 1, k1);
					  else //Boundary case
					     m_graph.add_arc(i0, i1, i2, k0, i0, i1, 0, k1 );
		
			
                } // for i2
            } // for i0
        } // for i1
    } // for i3

} // build_intersurface_arcs

////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////

} // namespace

#endif
