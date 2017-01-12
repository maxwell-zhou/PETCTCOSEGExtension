/* 
 ==========================================================================
 |   
 |   $Id: array_sorter.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___ARRAY_SORTER_HXX___
#   define ___ARRAY_SORTER_HXX___

#   include <optnet/config.h>
#   include <optnet/_base/array.hxx>
#   include <optnet/_base/array_ref.hxx>
#   include <algorithm>
#   include <vector>

/// @namespace optnet
namespace optnet {
    /// @namespace optnet::xtra
    namespace xtra {

//
// 3-D implementation
//
struct sort_index3
{
    size_t i[3];
};

// helper classes
        namespace detail {

            // helper: sort_asce3
            template <typename _Tx, typename _Tg>
            struct sort_asce3
            {
                sort_asce3(const optnet::array_base<_Tx, _Tg>& array)
                    : m_parray(&array)
                {
                }

                bool operator()(const sort_index2& i,
                                const sort_index2& j) const
                {
                    return (*m_parray)(i.i[0], i.i[1]) <=
                        (*m_parray)(j.i[0], j.i[1]);
                }

                optnet::array_base<_Tx, _Tg>* m_parray;
            };

            // helper: sort_desc3
            struct sort_desc3
            {
                sort_desc3(const optnet::array_base<_Tx, _Tg>& array)
                    : m_parray(&array)
                {
                }

                bool operator(const sort_index2& i,
                            const sort_index2& j) const
                {
                    return (*m_parray)(i.i[0], i.i[1]) >
                        (*m_parray)(j.i[0], j.i[1]);
                }

                optnet::array_base<_Tx, _Tg>* m_parray;
            };
        } // namespace

///////////////////////////////////////////////////////////////////////////
///  Sort the array indices in the specified order.
///////////////////////////////////////////////////////////////////////////
template <typename _Tx, typename _Tg>
void sort_array3(const optnet::array_base<_Tx, _Tg>& array, 
                 std::vector<sort_index3>&           vecindex,
                 bool                                descending = false
                 )
{
    using namespace optnet::utils::detail;
    
    if (descending) {
        std::sort(
            vecindex.begin(), 
            vecindex.end(), 
            sort_desc3<_Tx, _Tg>()
            );
    }
    else {
        std::sort(
            vecindex.begin(), 
            vecindex.end(), 
            sort_asce3<_Tx, _Tg>()
            );
    }
}

///////////////////////////////////////////////////////////////////////////
///  Sort the array indices in the specified order.
///////////////////////////////////////////////////////////////////////////
template <typename _Tx, typename _Tg>
void sort_array3(const optnet::array_base<_Tx, _Tg>& array, 
                 sort_index3*                        pindex,
                 size_t                              n,
                 bool                                descending = false
                 )
{
    using namespace optnet::utils::detail;
    
    if (descending) {
        std::sort(pindex, pindex + n, sort_desc3<_Tx, _Tg>());
    }
    else {
        std::sort(pindex, pindex + n, sort_asce3<_Tx, _Tg>());
    }
}

    }   // namespace
}   // namespace

#endif // ___ARRAY_SORTER_HXX___
