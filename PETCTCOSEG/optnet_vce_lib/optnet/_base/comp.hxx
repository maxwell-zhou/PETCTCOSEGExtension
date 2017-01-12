/*
 ==========================================================================
 |   
 |   $Id: comp.hxx 158 2005-02-09 23:37:29Z kangli $
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

#ifndef ___COMP_HXX___
#   define ___COMP_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <functional>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
template <typename _Ty>
struct roundoff_less : 
    public std::binary_function<_Ty, _Ty, bool>
{
    roundoff_less() :
        m_eps(_Ty()) {}

    roundoff_less(const _Ty& eps) :
        m_eps(eps) {}

    inline bool operator()(const _Ty& left,
                           const _Ty& right) const {
        return (left + m_eps < right);
    }

private:
    _Ty m_eps;
};

///////////////////////////////////////////////////////////////////////////
template <>
struct roundoff_less<double> : 
    public std::binary_function<double, double, bool>
{
    roundoff_less() :
        m_eps(1e-6) {}

    roundoff_less(const double& eps) :
        m_eps(eps) {}

    inline bool operator()(const double& left,
                           const double& right) const {
        return (left + m_eps < right);
    }

private:
    double m_eps;
};

///////////////////////////////////////////////////////////////////////////
template <>
struct roundoff_less<float> : 
    public std::binary_function<float, float, bool>
{
    roundoff_less() :
        m_eps(1e-6f) {}

    roundoff_less(const float& eps) :
        m_eps(eps) {}

    inline bool operator()(const float& left,
                           const float& right) const {
        return (left + m_eps < right);
    }

private:
    float m_eps;
};

} // namespace

#endif
