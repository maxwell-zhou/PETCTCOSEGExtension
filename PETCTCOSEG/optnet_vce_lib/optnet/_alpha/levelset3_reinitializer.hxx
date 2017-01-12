/*
 ==========================================================================
 |   
 |   $Id: levelset3_reinitializer.hxx 44 2005-01-19 03:18:25Z kangli $
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

#ifndef ___LEVELSET3_REINITIALIZER_HXX___
#   define ___LEVELSET3_REINITIALIZER_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4284)
#   endif

#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <optnet/_base/point3.hxx>
#   include <optnet/_alpha/levelset3_extractor.hxx>
#   include <optnet/_alpha/fast_marching_back3.hxx>
#   include <optnet/_alpha/fast_marching3.hxx>
#   include <cmath>

namespace optnet {

template <typename _Ty, typename _Tg = net_f_xy>
class levelset3_reinitializer
{
public:
    
    typedef _Ty                                     value_type;
    typedef value_type&                             reference;
    typedef const value_type&                       const_reference;

    typedef size_t                                  size_type;

    typedef fast_marching3<value_type,
                           double,
                           _Tg>                     fast_marching_type;
    typedef fast_marching_back3<value_type,
                                double,
                                _Tg>                fast_marching_back_type;

    typedef typename
        fast_marching_type::node_type               node_type;
    typedef typename
        fast_marching_type::node_vector_type        node_vector_type;

    typedef levelset3_extractor<value_type,
                                node_type,
                                node_vector_type,
                                _Tg>                levelset_extractor_type;

    typedef typename
        fast_marching_type::time_array_base_type    array_base_type;
    typedef typename
        fast_marching_type::time_array_ref_type     array_ref_type;
    typedef typename
        fast_marching_type::time_array_type         array_type;

    struct inner_indicator
    {
        inner_indicator(array_base_type&  levelset,
                        const value_type& levelset_value
                        ) :
            m_plevelset(&levelset),
            m_levelset_value(levelset_value)
        {
        }

        inline bool operator()(size_type i0,
                               size_type i1,
                               size_type i2
                               ) const
        {
            return ((*m_plevelset)(i0, i1, i2) >= m_levelset_value);
        }

    private:
        array_base_type*    m_plevelset;
        value_type          m_levelset_value;
    };

    struct outer_indicator
    {
        outer_indicator(array_base_type&  levelset,
                        const value_type& levelset_value
                        ) :
            m_plevelset(&levelset),
            m_levelset_value(levelset_value)
        {
        }

        inline bool operator()(size_type i0,
                               size_type i1,
                               size_type i2
                               ) const
        {
            return ((*m_plevelset)(i0, i1, i2) < m_levelset_value);
        }

    private:
        array_base_type*    m_plevelset;
        value_type          m_levelset_value;
    };


    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    ///////////////////////////////////////////////////////////////////////
    levelset3_reinitializer()
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Constructor.
    ///
    ///  @param levelset_value The value of the level set front.
    ///////////////////////////////////////////////////////////////////////
    levelset3_reinitializer(const value_type& levelset_value) :
        m_extractor(levelset_value)
    {
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Reinitialize the level set function to a signed distance function.
    ///
    ///  @param levelset The level set funtion.
    ///////////////////////////////////////////////////////////////////////
    void reinit(array_base_type& levelset)
    {
        node_vector_type& inner_nodes = m_fmm.trial_nodes();
        node_vector_type& outer_nodes = m_fmm_back.trial_nodes();

        m_extractor.extract(levelset, inner_nodes, outer_nodes);

        array_type levelset_copy(levelset);
        value_type levelset_value = get_levelset_value();

        m_fmm.create(levelset);
        m_fmm.solve(inner_indicator(levelset_copy, levelset_value));

        m_fmm_back.create(levelset);
        m_fmm_back.solve(outer_indicator(levelset_copy, levelset_value));
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Sets the value of the level set front.
    ///
    ///  @param levelset_value The value of the level set front.
    ///////////////////////////////////////////////////////////////////////
    inline void set_levelset_value(const value_type& levelset_value)
    {
        m_extractor.set_levelset_value(levelset_value);
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Returns the value of the level set front.
    ///
    ///  @return The value of the level set front.
    ///////////////////////////////////////////////////////////////////////
    inline const_reference get_levelset_value() const
    {
        return m_extractor.get_levelset_value();
    }


private:
    
    fast_marching_type      m_fmm;
    fast_marching_back_type m_fmm_back;
    levelset_extractor_type m_extractor;

};

} // namespace

#endif // ___LEVELSET3_REINITIALIZER_HXX___
