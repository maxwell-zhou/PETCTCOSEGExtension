/*
 ==========================================================================
 |   
 |   $Id: fast_marching_back3.hxx 186 2005-02-27 02:04:22Z kangli $
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

#ifndef ___FAST_MARCHING_BACK3_HXX___
#   define ___FAST_MARCHING_BACK3_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4127)
#   endif

#   include <optnet/_alpha/levelset_node.hxx>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/debug.hxx>
#   include <optnet/_base/point3.hxx>
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
///  @class fast_marching_back3
///  @brief 3-D backward fast marching method. This class marches the
///         arrival times backwards.
///////////////////////////////////////////////////////////////////////////
template <typename _Time,
          typename _Speed = double,
          typename _Tg = net_f_xy>
class fast_marching_back3
{
public:

    typedef _Time                                       time_value_type;
    typedef time_value_type&                            time_reference;
    typedef const time_value_type&                      time_const_reference;

    typedef _Speed                                      speed_value_type;
    typedef speed_value_type&                           speed_reference;
    typedef const speed_value_type&                     speed_const_reference;
    
    typedef unsigned char                               label_type;

    typedef array_base<time_value_type, _Tg>            time_array_base_type;
    typedef array_ref<time_value_type, _Tg>             time_array_ref_type;
    typedef array<time_value_type, _Tg>                 time_array_type;

    typedef array_base<time_value_type, _Tg>            speed_array_base_type;
    typedef array_ref<time_value_type, _Tg>             speed_array_ref_type;
    typedef array<time_value_type, _Tg>                 speed_array_type;

    typedef array_base<label_type, _Tg>                 label_array_base_type;
    typedef array_ref<label_type, _Tg>                  label_array_ref_type;
    typedef array<label_type, _Tg>                      label_array_type;

    typedef size_t                                      size_type;

    typedef point3<size_type>                           point_type;
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
    fast_marching_back3()
    {
        m_large_time = -static_cast<time_value_type>
            (std::numeric_limits<time_value_type>::max() / 2.0);
        m_speed_inv = -1.0;
    }

    ///////////////////////////////////////////////////////////////////////
    void create(time_array_base_type&        time_array,
                const speed_array_base_type* p_speed_array = NULL
                )
    {
        m_label_array.create_and_fill(FARAWAY,
                                      time_array.size_0(),
                                      time_array.size_1(),
                                      time_array.size_2()
                                      );
    
        // Save a pointer to the input time_array.
        m_p_time_array = &time_array;

        // Save a pointer of the speed array.
        m_p_speed_array = p_speed_array;
    }

    ///////////////////////////////////////////////////////////////////////
    inline node_vector_type& alive_nodes()
                { return m_alive_nodes; }
    inline const node_vector_type& alive_nodes() const
                { return m_alive_nodes; }

    ///////////////////////////////////////////////////////////////////////
    inline node_vector_type& trial_nodes()
                { return m_trial_nodes; }
    inline const node_vector_type& trial_nodes() const
                { return m_trial_nodes; }

    ///////////////////////////////////////////////////////////////////////
    void solve(double stopping_time = 1.0)
    {
        assert(NULL != m_p_time_array);

        time_value_type value;

        if (stopping_time > 0)
            stopping_time = m_large_time;

        m_p_time_array->fill(m_large_time);
        
        initialize_heap();
        
        while (!m_heap.empty()) {
            node_const_reference node = m_heap.top();
            m_heap.pop();

            value = (*m_p_time_array)(node.index);
            
            if (value != node.value) {
                debug::output
                    ("fast_marching_back3::solve: value != node.value\n");
                continue;
            }
            if (m_label_array(node.index) != TRIAL) {
                debug::output
                    ("fast_marching_back3::solve: m_label_array(node.index) != TRIAL\n",
                    m_label_array(node.index));
                continue;
            }
            if (value <= stopping_time) {
                debug::output
                    ("fast_marching_back3::solve: value < stopping_time\n");
                break;
            }

            m_label_array(node.index) = ALIVE;
            update_neighbors(node.index);

        } // while
    }

    ///////////////////////////////////////////////////////////////////////
    template <typename _RoiIndicator>
    void solve(_RoiIndicator roi,
               double        stopping_time = 1.0
               )
    {
        assert(NULL != m_p_time_array);

        time_value_type value;

        if (stopping_time > 0)
            stopping_time = m_large_time;

        for (size_type i2 = 0; i2 < m_p_time_array->size_2(); ++i2) {
            for (size_type i1 = 0; i1 < m_p_time_array->size_1(); ++i1) {
                for (size_type i0 = 0; i0 < m_p_time_array->size_0(); ++i0) {
                    if (roi(i0, i1, i2)) {
                        (*m_p_time_array)(i0, i1, i2) = m_large_time;
                    } // if
                }
            }
        }
        
        initialize_heap();
        
        while (!m_heap.empty()) {
            node_const_reference node = m_heap.top();
            m_heap.pop();

            value = (*m_p_time_array)(node.index);
            
            if (value != node.value) {
                debug::output
                    ("fast_marching_back3::solve: value != node.value\n");
                continue;
            }
            if (m_label_array(node.index) != TRIAL) {
                debug::output
                    ("fast_marching_back3::solve: m_label_array(node.index) != TRIAL\n",
                    m_label_array(node.index));
                continue;
            }
            if (value <= stopping_time) {
                debug::output
                    ("fast_marching_back3::solve: value < stopping_time\n");
                break;
            }

            m_label_array(node.index) = ALIVE;
            update_neighbors(roi, node.index);

        } // while
    }


private:

    typedef std::priority_queue<
        node_type,
        node_vector_type>               heap_type;

    ///////////////////////////////////////////////////////////////////////
    void initialize_heap()
    {
        node_const_iterator it_alive, it_trial;

        for (it_alive  = m_alive_nodes.begin();
             it_alive != m_alive_nodes.end(); ++it_alive) {

             node_const_reference node = *it_alive;
             size_type i0 = node.index[0];
             size_type i1 = node.index[1];
             size_type i2 = node.index[2];

             if (i0 > m_p_time_array->size_0() ||
                 i1 > m_p_time_array->size_1() ||
                 i2 > m_p_time_array->size_2()) {
                 // Ignore this node.
                 continue;
             }

             (*m_p_time_array)(i0, i1, i2) = node.value;
             m_label_array(i0, i1, i2) = ALIVE;
        } // for

        // Empty the trial heap if it is not.
        while (!m_heap.empty()) m_heap.pop();

        for (it_trial  = m_trial_nodes.begin();
             it_trial != m_trial_nodes.end(); ++it_trial) {

             node_const_reference node = *it_trial;
             size_type i0 = node.index[0];
             size_type i1 = node.index[1];
             size_type i2 = node.index[2];

             if (i0 > m_p_time_array->size_0() ||
                 i1 > m_p_time_array->size_1() ||
                 i2 > m_p_time_array->size_2()) {
                 // Ignore this node.
                 continue;
             }

             (*m_p_time_array)(i0, i1, i2) = node.value;
             m_label_array(i0, i1, i2) = TRIAL;
             m_heap.push(node);
        } // for
    }

    ///////////////////////////////////////////////////////////////////////
    bool update_neighbors(point_const_reference index)
    {
        size_type i, i0, i1, i2;

        i0 = index.v[0];
        i1 = index.v[1];
        i2 = index.v[2];
        
        //
        // six neighborhood
        //

        // i0-direction
        if (i0 > 0) {
            i = i0 - 1;
            if (m_label_array(i, i1, i2) != ALIVE) {
                update_node(i, i1, i2);
            }
        }
        if (i0 + 1 < m_label_array.size_0()) {
            i = i0 + 1;
            if (m_label_array(i, i1, i2) != ALIVE) {
                update_node(i, i1, i2);
            }
        }

        // i1-direction
        if (i1 > 0) {
            i = i1 - 1;
            if (m_label_array(i0, i, i2) != ALIVE) {
                update_node(i0, i, i2);
            }
        }
        if (i1 + 1 < m_label_array.size_1()) {
            i = i1 + 1;
            if (m_label_array(i0, i, i2) != ALIVE) {
                update_node(i0, i, i2);
            }
        }
        
        // i2-direction
        if (i2 > 0) {
            i = i2 - 1;
            if (m_label_array(i0, i1, i) != ALIVE) {
                update_node(i0, i1, i);
            }
        }
        if (i2 + 1 < m_label_array.size_2()) {
            i = i2 + 1;
            if (m_label_array(i0, i1, i) != ALIVE) {
                update_node(i0, i1, i);
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
        size_type i, i0, i1, i2;

        i0 = index.v[0];
        i1 = index.v[1];
        i2 = index.v[2];
        
        //
        // six neighborhood
        //

        // i0-direction
        if (i0 > 0) {
            i = i0 - 1;
            if (m_label_array(i, i1, i2) != ALIVE && roi(i, i1, i2)) {
                update_node(i, i1, i2);
            }
        }
        if (i0 + 1 < m_label_array.size_0()) {
            i = i0 + 1;
            if (m_label_array(i, i1, i2) != ALIVE && roi(i, i1, i2)) {
                update_node(i, i1, i2);
            }
        }

        // i1-direction
        if (i1 > 0) {
            i = i1 - 1;
            if (m_label_array(i0, i, i2) != ALIVE && roi(i0, i, i2)) {
                update_node(i0, i, i2);
            }
        }
        if (i1 + 1 < m_label_array.size_1()) {
            i = i1 + 1;
            if (m_label_array(i0, i, i2) != ALIVE && roi(i0, i, i2)) {
                update_node(i0, i, i2);
            }
        }
        
        // i2-direction
        if (i2 > 0) {
            i = i2 - 1;
            if (m_label_array(i0, i1, i) != ALIVE && roi(i0, i1, i)) {
                update_node(i0, i1, i);
            }
        }
        if (i2 + 1 < m_label_array.size_2()) {
            i = i2 + 1;
            if (m_label_array(i0, i1, i) != ALIVE && roi(i0, i1, i)) {
                update_node(i0, i1, i);
            }
        }

        return true;
    }

    ///////////////////////////////////////////////////////////////////////
    bool update_node(size_type i0, size_type i1, size_type i2)
    {
        static const size_type  _DIM_USED = 3;
        time_value_type         time[_DIM_USED];
        size_type               i;

        time[0] = m_large_time;
        time[1] = m_large_time;
        time[2] = m_large_time;

        // i0-direction
        if (i0 > 0) {
            i = i0 - 1;
            if (m_label_array(i, i1, i2) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i, i1, i2);
                time[0] = value;
            }
        }
        if (i0 + 1 < m_label_array.size_0()) {
            i = i0 + 1;
            if (m_label_array(i, i1, i2) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i, i1, i2);
                if (value > time[0]) {
                    time[0] = value;
                }
            }
        }

        // i1-direction
        if (i1 > 0) {
            i = i1 - 1;
            if (m_label_array(i0, i, i2) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i0, i, i2);
                time[1] = value;
            }
        }
        if (i1 + 1 < m_label_array.size_1()) {
            i = i1 + 1;
            if (m_label_array(i0, i, i2) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i0, i, i2);
                if (value > time[1]) {
                    time[1] = value;
                }
            }
        }

        // i2-direction
        if (i2 > 0) {
            i = i2 - 1;
            if (m_label_array(i0, i1, i) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i0, i1, i);
                time[2] = value;
            }
        }
        if (i2 + 1 < m_label_array.size_2()) {
            i = i2 + 1;
            if (m_label_array(i0, i1, i) == ALIVE) {
                time_const_reference value 
                    = (*m_p_time_array)(i0, i1, i);
                if (value > time[2]) {
                    time[2] = value;
                }
            }
        }

        //
        // Algorithm: see fast_marching.hxx
        //

        double b, c, solution, discrim;
        int    a = 0;

        // Solve the quadratic equation: Equation 8 in Reference [1].
        solution = (double)m_large_time;
        b        = 0.0;
        
        if (NULL != m_p_speed_array) {
            c = (double)(*m_p_speed_array)(i0, i1, i2);
            c = -1.0 / (c * c);
        }
        else {
            c = m_speed_inv;
        }

        // Sort the local node list (in nonascending order).
        std::sort(&(time[0]), &(time[0]) + _DIM_USED, 
            std::greater<time_value_type>());

        for (i = 0; i < _DIM_USED; ++i) {

           if (solution > time[i]) break;

            a += 1;
            b += time[i];
            c += time[i] * time[i];

            discrim = b * b - a * c;
            
            if (discrim < 0.0) {
                debug::output(
                    "fast_marching_back3::update_node: No real-valued solution."
                    );
                return false;
            }

            solution = (-sqrt(discrim) + b) / a;
        }


        if (solution > (double)m_large_time) {
            node_type node;

            assert(m_label_array(i0, i1, i2) != ALIVE);

            time_value_type value
                = static_cast<time_value_type>(solution);

            if (value > (*m_p_time_array)(i0, i1, i2))
                (*m_p_time_array)(i0, i1, i2) = value;

            // Insert the point into the trail heap.
            m_label_array(i0, i1, i2) = TRIAL;

            node.value    = value;
            node.index[0] = i0;
            node.index[1] = i1;
            node.index[2] = i2;

            m_heap.push(node);
        }

        return true;
    }


    ///////////////////////////////////////////////////////////////////////
    time_value_type                 m_large_time;
    double                          m_speed_inv;
    
    label_array_type                m_label_array;
    const speed_array_base_type*    m_p_speed_array;
    time_array_base_type*           m_p_time_array;

    // node containers
    node_vector_type                m_trial_nodes;
    node_vector_type                m_alive_nodes;
    heap_type                       m_heap;
};

// constants
template <typename _Time, typename _Speed, typename _Tg>
const typename fast_marching_back3<_Time, _Speed, _Tg>::label_type
    fast_marching_back3<_Time, _Speed, _Tg>::FARAWAY = 
        (typename fast_marching_back3<_Time, _Speed, _Tg>::label_type)(0);

template <typename _Time, typename _Speed, typename _Tg>
const typename fast_marching_back3<_Time, _Speed, _Tg>::label_type
    fast_marching_back3<_Time, _Speed, _Tg>::ALIVE = 
        (typename fast_marching_back3<_Time, _Speed, _Tg>::label_type)(1);

template <typename _Time, typename _Speed, typename _Tg>
const typename fast_marching_back3<_Time, _Speed, _Tg>::label_type
    fast_marching_back3<_Time, _Speed, _Tg>::TRIAL = 
        (typename fast_marching_back3<_Time, _Speed, _Tg>::label_type)(2);

} // namespace

#endif // ___FAST_MARCHING_BACK3_HXX___
