/*
==========================================================================
|   
|   $Id: iterator.hxx 186 2005-02-27 02:04:22Z kangli $
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

#ifndef ___ITERATOR_HXX___
#   define ___ITERATOR_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#   endif

#   include <iterator>
#   include <optnet/config.h>

/// @namespace optnet
namespace optnet {

namespace detail {

//
// class normal_iterator
//
#   ifndef __OPTNET_CRAPPY_MSC__

    ///////////////////////////////////////////////////////////////////////
    ///  @class normal_iterator
    ///  @brief Standard iterator type.
    ///////////////////////////////////////////////////////////////////////
    template<typename _Iterator, typename _Container>
    class normal_iterator
        : public std::iterator<
            typename std::iterator_traits<_Iterator>::iterator_category,
            typename std::iterator_traits<_Iterator>::value_type,
            typename std::iterator_traits<_Iterator>::difference_type,
            typename std::iterator_traits<_Iterator>::pointer,
            typename std::iterator_traits<_Iterator>::reference>
    {
    protected:
        _Iterator _M_current;

    public:
        typedef typename std::iterator_traits<_Iterator>::difference_type    
            difference_type;
        typedef typename std::iterator_traits<_Iterator>::reference reference;
        typedef typename std::iterator_traits<_Iterator>::pointer   pointer;

        normal_iterator() : _M_current(_Iterator()) { }

        explicit 
            normal_iterator(const _Iterator& __i) : _M_current(__i) { }

        // Allow iterator to const_iterator conversion
        template<typename _Iter>
        inline normal_iterator(const normal_iterator<_Iter, _Container>& __i)
        : _M_current(__i.base()) { }

        // Forward iterator requirements
        reference operator* () const { return *_M_current; }
        pointer   operator->() const { return  _M_current; }

        normal_iterator& operator++() { ++_M_current; return *this; }
        normal_iterator  operator++(int)
        { return normal_iterator(_M_current++); }

        // Bidirectional iterator requirements
        normal_iterator& operator--() { --_M_current; return *this; }
        normal_iterator  operator--(int)
        { return normal_iterator(_M_current--); }

        // Random access iterator requirements
        reference        operator[](const difference_type& __n) const
        { return _M_current[__n]; }

        normal_iterator& operator+=(const difference_type& __n)
        { _M_current += __n; return *this; }

        normal_iterator  operator+ (const difference_type& __n) const
        { return normal_iterator(_M_current + __n); }

        normal_iterator& operator-=(const difference_type& __n)
        { _M_current -= __n; return *this; }

        normal_iterator  operator- (const difference_type& __n) const
        { return normal_iterator(_M_current - __n); }

        const _Iterator& base() const { return _M_current; }
    };

    // Forward iterator requirements
    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator==(const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() == __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator==(const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() == __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator!=(const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() != __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator!=(const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() != __rhs.base(); }

    // Random access iterator requirements
    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool 
    operator< (const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() < __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator< (const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() < __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator> (const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() > __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator> (const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() > __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator<=(const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() <= __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator<=(const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() <= __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline bool
    operator>=(const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() >= __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline bool
    operator>=(const normal_iterator<_Iterator,  _Container>& __lhs,
               const normal_iterator<_Iterator,  _Container>& __rhs)
    { return __lhs.base() >= __rhs.base(); }

    template<typename _IteratorL, typename _IteratorR, typename _Container>
    inline typename normal_iterator<_IteratorL, _Container>::difference_type
    operator- (const normal_iterator<_IteratorL, _Container>& __lhs,
               const normal_iterator<_IteratorR, _Container>& __rhs)
    { return __lhs.base() - __rhs.base(); }

    template<typename _Iterator, typename _Container>
    inline normal_iterator<_Iterator, _Container>
    operator+ (
        typename normal_iterator<_Iterator, _Container>::difference_type __n,
        const    normal_iterator<_Iterator, _Container>& __i)
    { return normal_iterator<_Iterator, _Container>(__i.base() + __n); }

#   endif // __OPTNET_CRAPPY_MSC__

} // namespace

} // namespace


#endif 
