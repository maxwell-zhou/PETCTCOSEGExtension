/*
 ==========================================================================
 |   
 |   $Id: fast_marching2.hxx 186 2005-02-27 02:04:22Z kangli $
 |
 |   Written by Kang Li <kangl@cmu.edu>
 |   Department of Electrical and Computer Engineering
 |   Carnegie Mellon University
 |   
 ==========================================================================
 |   This file is a part of the OptimalNet library.
 ==========================================================================
 | Copyright (c) 2004-2005 Kang Li <kangl@cmu.edu>. All Rights Reserved.
 | 
 | This software is supplied under the terms of a license agreement or
 | nondisclosure agreement  with the author  and may not be copied  or
 | disclosed except in accordance with the terms of that agreement.
 ==========================================================================
 */

#ifndef ___FAST_MARCHING2_HXX___
#   define ___FAST_MARCHING2_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4127)
#   endif

#   include <optnet/_alpha/levelset_node.hxx>
#   include <optnet/_base/array2.hxx>
#   include <optnet/_base/array2_ref.hxx>
#   include <optnet/_base/debug.hxx>
#   include <optnet/_base/point2.hxx>
#   include <functional>
#   include <limits>
#   include <queue>
#   include <vector>

#   ifdef max
#       undef max
#   endif

/// @namespace optnet
namespace optnet {

///////////////////////////////////////////////////////////////////////////
///  @class fast_marching2
///  @brief 2-D fast marching method.
///////////////////////////////////////////////////////////////////////////
template <typename _Time, typename _Speed = double>
class fast_marching2
{
public:

    typedef _Time                                       time_value_type;
    typedef time_value_type&                            time_reference;
    typedef const time_value_type&                      time_const_reference;

    typedef _Speed                                      speed_value_type;
    typedef speed_value_type&                           speed_reference;
    typedef const speed_value_type&                     speed_const_reference;
    
    typedef unsigned char                               label_type;

    typedef array2_base<time_value_type>                time_array_base_type;
    typedef array2_ref<time_value_type>                 time_array_ref_type;
    typedef array2<time_value_type>                     time_array_type;

    typedef array2_base<time_value_type>                speed_array_base_type;
    typedef array2_ref<time_value_type>                 speed_array_ref_type;
    typedef array2<time_value_type>                     speed_array_type;

    typedef array2_base<label_type>                     label_array_base_type;
    typedef array2_ref<label_type>                      label_array_ref_type;
    typedef array2<label_type>                          label_array_type;

    typedef size_t                                      size_type;

    typedef point2<size_type>                           point_type;
    typedef point_type&                                 point_reference;
    typedef const point_type&                           point_const_reference;

    typedef levelset_node<time_value_type,
                          point_type>                   node_type;

    typedef std::vector<node_type>                      node_vector_type;
    typedef typename node_vector_type::reference        node_reference;
    typedef typename node_vector_type::const_reference  node_const_reference;
    typedef typename node_vector_type::iterator         node_iterator;
    typedef typename node_vector_type::const_iterator   node_const_iterator;
    typedef typename node_vector_type::pointer          node_pointer;
    typedef typename node_vector_type::const_pointer    node_const_pointer;

    static const label_type FARAWAY;    /// Faraway node.
    static const label_type ALIVE;      /// Alive node.
    static const label_type TRIAL;      /// Trial node.

    ///////////////////////////////////////////////////////////////////////
    /// Constructor.
    ///////////////////////////////////////////////////////////////////////
    fast_marching2()
    {
        m_large_time = static_cast<time_value_type>
            (std::numeric_limits<time_value_type>::max() / 2.0);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Initialize the fast marching solver.
    ///
    ///  @param[out] time_array     The output time array containing the
    ///                             arrival time at each point.
    ///  @param[in]  p_speed_array  The input speed array.
    ///
    ///////////////////////////////////////////////////////////////////////
    void create(time_array_base_type&        time_array,
                const speed_array_base_type* p_speed_array = NULL
                )
    {
        m_label_array.create_and_fill(FARAWAY,
                                      time_array.size_0(),
                                      time_array.size_1()
                                      );
    
        // Save a pointer of the speed array.
        m_p_speed_array = p_speed_array;

        // Save a pointer to the time_array.
        m_p_time_array = &time_array;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a reference to the vector of ALIVE nodes.
    ///////////////////////////////////////////////////////////////////////
    inline node_vector_type& alive_nodes()
                { return m_alive_nodes; }
    ///////////////////////////////////////////////////////////////////////
    ///  Returns a constant reference to the vector of ALIVE nodes.
    ///////////////////////////////////////////////////////////////////////
    inline const node_vector_type& alive_nodes() const
                { return m_alive_nodes; }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns a reference to the vector of TRIAL nodes.
    ///////////////////////////////////////////////////////////////////////
    inline node_vector_type& trial_nodes()
                { return m_trial_nodes; }
    ///////////////////////////////////////////////////////////////////////
    ///  Returns a constant reference to the vector of TRIAL nodes.
    ///////////////////////////////////////////////////////////////////////
    inline const node_vector_type& trial_nodes() const
                { return m_trial_nodes; }

    ///////////////////////////////////////////////////////////////////////
    ///  Starts fast marching.
    ///
    ///  @param stopping_time  Stop marching when the arrival time reaches
    ///                        this value.
    ///////////////////////////////////////////////////////////////////////
    void solve(double stopping_time = -1.0)
    {
        assert(NULL != m_p_time_array);

        time_value_type value;

        if (stopping_time < 0)
            stopping_time = m_large_time;

        m_p_time_array->fill(m_large_time);

        initialize_heap();
        
        while (!m_heap.empty()) {
            node_const_reference node = m_heap.top();
            m_heap.pop();

            value = (*m_p_time_array)(node.index);
            
            if (value != node.value) {
                debug::output
                    ("fast_marching2::solve: value != node.value\n");
                continue;
            }
            if (m_label_array(node.index) != TRIAL) {
                debug::output
                    ("fast_marching2::solve: m_label_array(node.index) != TRIAL\n",
                    m_label_array(node.index));
                continue;
            }
            if (value >= stopping_time) {
                debug::output
                    ("fast_marching2::solve: value > stopping_time\n");
                break;
            }

            m_label_array(node.index) = ALIVE;
            update_neighbors(node.index);

        } // while
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Starts fast marching in the specified region-of-interest (ROI).
    ///
    ///  @param roi            The region-of-interest indicator functor,
    ///                        where roi(x, y) returns true if the point
    ///                        (x, y) is in the ROI.
    ///  @param stopping_time  Stop marching when the arrival time reaches
    ///                        this value.
    ///////////////////////////////////////////////////////////////////////
    template <typename _RoiIndicator>
    void solve(_RoiIndicator roi,
               double        stopping_time = -1.0
               )
    {
        assert(NULL != m_p_time_array);

        time_value_type value;

        if (stopping_time < 0)
            stopping_time = m_large_time;

        for (size_type i1 = 0; i1 < m_p_time_array->size_1(); ++i1) {
            for (size_type i0 = 0; i0 < m_p_time_array->size_0(); ++i0) {
                if (roi(i0, i1)) {
                    (*m_p_time_array)(i0, i1) = m_large_time;
                } // if
            }
        }

        initialize_heap();
        
        while (!m_heap.empty()) {
            node_const_reference node = m_heap.top();
            m_heap.pop();

            value = (*m_p_time_array)(node.index);
            
            if (value != node.value) {
                debug::output
                    ("fast_marching2::solve: value != node.value\n");
                continue;
            }
            if (m_label_array(node.index) != TRIAL) {
                debug::output
                    ("fast_marching2::solve: m_label_array(node.index) != TRIAL\n",
                    m_label_array(node.index));
                continue;
            }
            if (value > stopping_time) {
                debug::output
                    ("fast_marching2::solve: value > stopping_time\n");
                break;
            }

            m_label_array(node.index) = ALIVE;
            update_neighbors(roi, node.index);

        } // while
    }


private:

    typedef std::greater<node_type>     node_comparer;
    typedef std::priority_queue<
        node_type,
        node_vector_type,
        node_comparer>                  heap_type;

    ///////////////////////////////////////////////////////////////////////
    void initialize_heap()
    {
        node_const_iterator it_alive, it_trial;

        for (it_alive  = m_alive_nodes.begin();
             it_alive != m_alive_nodes.end(); ++it_alive) {

             node_const_reference node = *it_alive;
             size_type i0 = node.index[0];
             size_type i1 = node.index[1];

             if (i0 > m_p_time_array->size_0() ||
                 i1 > m_p_time_array->size_1()) {
                 // Ignore this node.
                 continue;
             }

             (*m_p_time_array)(i0, i1) = node.value;
             m_label_array(i0, i1) = ALIVE;
        } // for

        // Empty the trial heap if it is not.
        while (!m_heap.empty()) m_heap.pop();

        for (it_trial  = m_trial_nodes.begin();
             it_trial != m_trial_nodes.end(); ++it_trial) {

             node_const_reference node = *it_trial;
             size_type i0 = node.index[0];
             size_type i1 = node.index[1];

             if (i0 > m_p_time_array->size_0() ||
                 i1 > m_p_time_array->size_1()) {
                 // Ignore this node.
                 continue;
             }

             (*m_p_time_array)(i0, i1) = node.value;
             m_label_array(i0, i1) = TRIAL;
             m_heap.push(node);
        } // for
    }

    ///////////////////////////////////////////////////////////////////////
    bool update_neighbors(point_const_reference index)
    {
        size_type i, i0, i1;

        i0 = index.v[0];
        i1 = index.v[1];
        
        // i0-direction
        if (i0 > 0) {
            i = i0 - 1;
            if (m_label_array(i, i1) != ALIVE) {
                update_node(i, i1);
            }
        }
        if (i0 + 1 < m_label_array.size_0()) {
            i = i0 + 1;
            if (m_label_array(i, i1) != ALIVE) {
                update_node(i, i1);
            }
        }

        // i1-direction
        if (i1 > 0) {
            i = i1 - 1;
            if (m_label_array(i0, i) != ALIVE) {
                update_node(i0, i);
            }
        }
        if (i1 + 1 < m_label_array.size_1()) {
            i = i1 + 1;
            if (m_label_array(i0, i) != ALIVE) {
                update_node(i0, i);
            }
        }
        
        return true;
    }

    ///////////////////////////////////////////////////////////////////////
    template <typename _RoiIndicator>
    bool update_neighbors(_RoiIndicator         roi,
                          point_const_reference index
                          )
    {
        size_type i, i0, i1;

        i0 = index.v[0];
        i1 = index.v[1];
        
        // i0-direction
        if (i0 > 0) {
            i = i0 - 1;
            if (m_label_array(i, i1) != ALIVE && roi(i, i1)) {
                update_node(i, i1);
            }
        }
        if (i0 + 1 < m_label_array.size_0()) {
            i = i0 + 1;
            if (m_label_array(i, i1) != ALIVE && roi(i, i1)) {
                update_node(i, i1);
            }
        }

        // i1-direction
        if (i1 > 0) {
            i = i1 - 1;
            if (m_label_array(i0, i) != ALIVE && roi(i0, i)) {
                update_node(i0, i);
            }
        }
        if (i1 + 1 < m_label_array.size_1()) {
            i = i1 + 1;
            if (m_label_array(i0, i) != ALIVE && roi(i0, i)) {
                update_node(i0, i);
            }
        }
        
        return true;
    }

    ///////////////////////////////////////////////////////////////////////
    bool update_node(size_type i0, size_type i1)
    {
        size_type       i;
        node_type       node;
        time_value_type time[2];

        time[0] = m_large_time;
        time[1] = m_large_time;

        // i0-direction
        if (i0 > 0) {
            i = i0 - 1;
            if (m_label_array(i, i1) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i, i1);
                time[0] = value;
            }
        }
        if (i0 + 1 < m_label_array.size_0()) {
            i = i0 + 1;
            if (m_label_array(i, i1) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i, i1);
                if (value < time[0]) {
                    time[0] = value;
                }
            }
        }

        // i1-direction
        if (i1 > 0) {
            i = i1 - 1;
            if (m_label_array(i0, i) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i0, i);
                time[1] = value;
            }
        }
        if (i1 + 1 < m_label_array.size_1()) {
            i = i1 + 1;
            if (m_label_array(i0, i) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i0, i);
                if (value < time[1]) {
                    time[1] = value;
                }
            }
        }

        double b, p, solution, discrim;

        // Solve the quadratic equation: Equation 8 in Reference [1].
        solution = (double)m_large_time;
        
        if (NULL != m_p_speed_array) {
            p = 1.0 / (double)(*m_p_speed_array)(i0, i1);
        }
        else {
            p = 1.0;
        }

        if (time[0] > time[1]) {
            time_value_type tmp = time[0];
            time[0] = time[1];
            time[1] = tmp;
        }

        b = time[1] - time[0];
        b = b * b;

        if (p * p > b) {
            
            discrim = 2 * p * p - b;

            if (discrim < 0.0) {
                debug::output(
                    "fast_marching2::update_node: No real-valued solution."
                    );
                return false;
            }
            
            solution = (sqrt(discrim) + time[0] + time[1]) / 2.0;
        }
        else
            solution = time[0] + p;

        if (solution < (double)m_large_time) {
            assert(m_label_array(i0, i1) != ALIVE);

            time_value_type value
                = static_cast<time_value_type>(solution);

            if (value < (*m_p_time_array)(i0, i1))
                (*m_p_time_array)(i0, i1) = value;

            // Insert the point into the trail heap.
            m_label_array(i0, i1) = TRIAL;

            node.value    = value;
            node.index[0] = i0;
            node.index[1] = i1;

            m_heap.push(node);
        }

        return true;
    }


    ///////////////////////////////////////////////////////////////////////
    time_value_type                 m_large_time;
    
    label_array_type                m_label_array;
    const speed_array_base_type*    m_p_speed_array;
    time_array_base_type*           m_p_time_array;

    // node containers
    node_vector_type                m_trial_nodes;
    node_vector_type                m_alive_nodes;
    heap_type                       m_heap;
};

// constants
template <typename _Time, typename _Speed>
const typename fast_marching2<_Time, _Speed>::label_type
    fast_marching2<_Time, _Speed>::FARAWAY = 
        (typename fast_marching2<_Time, _Speed>::label_type)(0);

template <typename _Time, typename _Speed>
const typename fast_marching2<_Time, _Speed>::label_type
    fast_marching2<_Time, _Speed>::ALIVE = 
        (typename fast_marching2<_Time, _Speed>::label_type)(1);

template <typename _Time, typename _Speed>
const typename fast_marching2<_Time, _Speed>::label_type
    fast_marching2<_Time, _Speed>::TRIAL = 
        (typename fast_marching2<_Time, _Speed>::label_type)(2);

} // namespace

#endif // ___FAST_MARCHING2_HXX___
