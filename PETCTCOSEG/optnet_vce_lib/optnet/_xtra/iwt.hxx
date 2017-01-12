/* 
 ==========================================================================
 |   
 |   $Id: iwt.hxx 186 2005-02-27 02:04:22Z kangli $
 |
 |   Written by Kang Li <likang@ieee.org>
 |   Department of Electrical and Computer Engineering
 |   Carnegie Mellon University
 |   
 ==========================================================================
 |   This file is a part of the OptimalNet library.
 ==========================================================================
 | Copyright (c) 2004-2005 Kang Li <likang@ieee.org>. All Rights Reserved.
 | 
 | This software is supplied under the terms of a license agreement or
 | nondisclosure agreement  with the author  and may not be copied  or
 | disclosed except in accordance with the terms of that agreement.
 ==========================================================================
 */

#ifndef ___IWT_HXX___
#   define ___IWT_HXX___

/*
 ==========================================================================
  - Purpose:

      This file implements the Interactive Watershed Transform (IWT).

  - References:
    [1] Horst K. Hahn and Heinz-Otto Peitgen
        IWT ¨C Interactive Watershed Transform: A hierarchical method for
          efficient interactive and automated segmentation of
          multidimensional grayscale images,
        Proc. Medical Imaging, SPIE 5032, San Diego, Feb 2003.
    [2] VINCENT, L., and SOILLE, P. 1991.
        Watersheds in digital spaces: An efficient algorithm based on
          immersion simulations.
        IEEE Transactions on Pattern Analysis and Machine Intelligence
        Vol. 13, Iss. 6 (June), 583?98.
    [3] Jos B.T.M Boerdink and Arnold Meijster
        The watershed transform: Definitions, algorithms and
          parallelization strategies
        Fundamenta Informaticae, 41 (2001) 187-228.
 ==========================================================================
 */

#   include <optnet/_xtra/array_sorter.hxx>

/// @namespace optnet
namespace optnet {
    /// @namespace optnet::xtra
    namespace xtra {

///////////////////////////////////////////////////////////////////////////
/// @class iwt
/// @brief A class that performs the interactive watershed transform.
///////////////////////////////////////////////////////////////////////////
template <typename _Tx,
          typename _Tg = net_f_xy>
class iwt
{
    typedef optnet::array_base<_Tx, _Tg>    array_base_type;
    typedef optnet::array_ref <_Tx, _Tg>    array_ref_type;
    typedef optnet::array     <_Tx, _Tg>    array_type;

    typedef _Tx                             value_type; 
    typedef size_t                          size_type;

public:

    iwt(const array_base_type& rarray)
    {
        //TODO: Implement this.
    }

    void compute()
    {
        //TODO: Implement this.
    }
    
};

    }   // namespace
}   // namespace

#endif // ___IWT_HXX___
