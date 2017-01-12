/*
 ==========================================================================
 |   
 |   $Id: graphspec.hxx 177 2005-02-17 20:26:01Z kangli $
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

#ifndef ___GRAPHSPEC_HXX___
#   define ___GRAPHSPEC_HXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4284)
#       pragma warning(disable: 4786)
#       pragma warning(disable: 4571)
#   endif

#   include <optnet/config.h>
#   include <optnet/_base/point.hxx>
#   include <optnet/_utils/xstring.hxx>
#   include <algorithm>
#   include <fstream>
#   include <vector>
#   include <cstring>

/// @namespace optnet
namespace optnet {

#   define OPTNET_GRAPHSPEC_MAX_SURF 2

///////////////////////////////////////////////////////////////////////////
///  @class graphspec
///  @brief Graph specification.
///////////////////////////////////////////////////////////////////////////
class graphspec
{
public:
    
    /// @class column_type
    /// @brief Specifies a column in a graph
    struct column_type
    {
        int     id;                                     // unique id
        double  pos[3 * OPTNET_GRAPHSPEC_MAX_SURF];     // position
        double  normal[3 * OPTNET_GRAPHSPEC_MAX_SURF];  // normal vector
        double  fact[OPTNET_GRAPHSPEC_MAX_SURF];        // scaling factor

        column_type& operator+=(const column_type& rhs)
        {
            int i;
            for (i = 0; i < sizeof(pos) / sizeof(double); ++i)
                pos[i] += rhs.pos[i];
            for (i = 0; i < sizeof(normal) / sizeof(double); ++i)
                normal[i] += rhs.normal[i];
            return *this;
        }
        
        column_type& operator/=(const double& rhs)
        {
            int i;
            for (i = 0; i < sizeof(pos) / sizeof(double); ++i)
                pos[i] /= rhs;
            for (i = 0; i < sizeof(normal) / sizeof(double); ++i)
                normal[i] /= rhs;
            return *this;
        }
    };

    /// @class adjacency_type
    /// @brief Defines column adjacencies in a graph
    struct adjacency_type
    {
        int     id[2];

        bool operator<(const adjacency_type& rhs) const
            { return (memcmp(id, rhs.id, sizeof(int) * 2) < 0); }
    };

    /// @class triplet_type
    /// @brief Defines column triplets in a graph
    struct triplet_type
    {
        int     id[3];
    };

    // column comparator
    struct column_cmp
    {
        bool operator()(const column_type& lhs,
                        const column_type& rhs) {
            return (lhs.id < rhs.id);
        }
    };

    typedef column_type&            column_reference;
    typedef const column_type&      column_const_reference;
    typedef adjacency_type&         adjacency_reference;
    typedef const adjacency_type&   adjacency_const_reference;
    typedef triplet_type&           triplet_reference;
    typedef const triplet_type&     triplet_const_reference;

    typedef point<double, 3>        point_type;
    typedef point_type&             point_reference;
    typedef const point_type&       point_const_reference;

    typedef point<double, 6>        boundary_type;
    typedef boundary_type&          boundary_reference;
    typedef const boundary_type&    boundary_const_reference;

    typedef size_t                  size_type;


    ///////////////////////////////////////////////////////////////////////
    ///  Default constructor.
    graphspec() : 
        m_sorted(false),
        m_column_type(1),
        m_adjacency_type(1),
        m_triplet_type(1)
    {
    }
    
    ///////////////////////////////////////////////////////////////////////
    /// Set sorted flag.
    inline void set_sorted(bool sorted) { m_sorted = sorted; }

    ///////////////////////////////////////////////////////////////////////
    ///  Sort columns by id in ascending order.
    inline void sort_columns()
    {
        std::sort(m_columns.begin(), m_columns.end(), column_cmp()); 
        m_sorted = true; // set sorted flag
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Gets the column type.
    inline int  get_column_type() const      { return m_column_type; }
    ///////////////////////////////////////////////////////////////////////
    ///  Sets the column type.
    inline void set_column_type(int type)    { m_column_type = type; }

    ///////////////////////////////////////////////////////////////////////
    ///  The number of columns in the graph.
    inline size_type num_columns() const     { return m_columns.size(); }
    ///  The number of adjacency specifications.
    inline size_type num_adjacencies() const { return m_adjacencies.size(); }
    ///  The number of triplet specifications.
    inline size_type num_triplets() const    { return m_triplets.size(); }

    ///////////////////////////////////////////////////////////////////////
    ///  Get column element.
    inline column_const_reference get_column(size_type i) const
                    { return m_columns[i]; }
    ///  Get adjacency element.
    inline adjacency_const_reference get_adjacency(size_type i) const
                    { return m_adjacencies[i]; }
    ///  Get triplet element.
    inline triplet_const_reference get_triplet(size_type i) const
                    { return m_triplets[i]; }
    ///  Get boundary.
    inline boundary_const_reference get_boundary() const
                    { return m_boundary; }
    ///  Get origin.
    inline point_const_reference get_origin() const
                    { return m_origin; }

    ///////////////////////////////////////////////////////////////////////
    ///  Set column element.
    inline void set_column(size_type i, column_const_reference c)
    {
        m_columns[i] = c;
        if (c.id != (int)i) m_sorted = false;
    }
    ///  Set adjacency element.
    inline void set_adjacency(size_type i, adjacency_const_reference a)
                    { m_adjacencies[i] = a; }
    ///  Set triplet element.
    inline void set_triplet(size_type i, triplet_const_reference t)
                    { m_triplets[i] = t; }
    ///  Set boundary.
    inline void set_boundary(boundary_const_reference b) { m_boundary = b; }
    inline void set_boundary(double bleft,
                             double bright,
                             double bbottom,
                             double btop,
                             double bnear,
                             double bfar)
    {
        m_boundary[0] = bleft;
        m_boundary[1] = bright;
        m_boundary[2] = bbottom;
        m_boundary[3] = btop;
        m_boundary[4] = bnear;
        m_boundary[5] = bfar;
    }

    ///  Set origin.
    inline void set_origin(point_const_reference o) { m_origin = o; }
    inline void set_origin(double x, double y, double z)
    {
        m_origin[0] = x;
        m_origin[1] = y;
        m_origin[2] = z;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Add a graph column.
    inline void add_column(column_const_reference c)
    {
        if (c.id != (int)(&(m_columns.back()) - &(m_columns.front())))
            m_sorted = false;
        m_columns.push_back(c);
    }
    ///  Add an adjacency.
    inline void add_adjacency(adjacency_const_reference a)
                    { m_adjacencies.push_back(a); }
    ///  Add a triplet.
    inline void add_triplet(triplet_const_reference t)
                    { m_triplets.push_back(t); }

    // convenient functions

    /// Get scaling factor for column i.
    inline double get_fact(size_type i, size_type j) const
    {
        return m_columns[i].fact[j];
    }

    /// Set scaling factor for column i.
    inline void set_fact(size_type i, size_type j, const double& fact)
    {
        m_columns[i].fact[j] = fact;
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Create adjacencies using information of the triplets.
    ///////////////////////////////////////////////////////////////////////
    inline void adjacencies_from_triplets();

    ///////////////////////////////////////////////////////////////////////
    ///  Recalculate column positions from scaling factors.
    ///////////////////////////////////////////////////////////////////////
    inline void positions_from_facts();

    ///////////////////////////////////////////////////////////////////////
    ///  Make the target mesh smoother.
    ///////////////////////////////////////////////////////////////////////
    inline void smooth_mesh();

    ///////////////////////////////////////////////////////////////////////
    ///  Clear all loaded data.
    inline void clear()
    {
             m_origin.clear();
           m_boundary.clear();
            m_columns.clear();
        m_adjacencies.clear();
           m_triplets.clear();
            m_sorted = false;
    }

    ///////////////////////////////////////////////////////////////////////
    inline void load(const char* filename)
    {
        using namespace optnet::utils;

        if (str_ends_with_no_case(filename, ".ons")) {
            load_source(filename);
        }
        else if (str_ends_with_no_case(filename, ".ont")) {
            load_target(filename);
        }
        else {
            bool error_occurred = false;

            try {
                load_source(filename);
            } catch (...) {
                error_occurred = true;
            }
            if (!error_occurred)
                return;

            try {
                load_target(filename);
            } catch (...) {
                error_occurred = true;
            }
            if (!error_occurred)
                return;

            char errmsg[4096];
            secure_sprintf(errmsg, 4096,
                "graphspec::load: Could not load %s.", filename);
            throw_exception(std::runtime_error(errmsg));
        }
    }

    ///////////////////////////////////////////////////////////////////////
    inline void save(const char* filename)
    {
        using namespace optnet::utils;

        if (str_ends_with_no_case(filename, ".ons")) {
            save_source(filename);
        }
        else if (str_ends_with_no_case(filename, ".ont")) {
            save_target(filename);
        }
        else {
            char errmsg[4096];
            secure_sprintf(errmsg, 4096, "graphspec::save:"
                "Please specify an extension for the file name '%s'.",
                filename);
            throw_exception(std::runtime_error(errmsg));
        }
    }

    ///////////////////////////////////////////////////////////////////////
    ///  Load the graph source specification from an xml file.
    inline void load_source(const char* filename);
    ///  Save the graph source specification to an xml file.
    inline void save_source(const char* filename);
    ///  Save the graph source specification to a stream.
    inline void save_source(std::ostream& os);

    ///////////////////////////////////////////////////////////////////////
    ///  Load the graph target specification from an xml file.
    inline void load_target(const char* filename);
    ///  Save the graph target specification to an xml file.
    inline void save_target(const char* filename);
    ///  Save the graph target specification to a stream.
    inline void save_target(std::ostream& os);


private:
    
    typedef std::vector<column_type>    column_vector;
    typedef std::vector<adjacency_type> adjacency_vector;
    typedef std::vector<triplet_type>   triplet_vector;

    bool                m_sorted;
    int                 m_column_type;
    int                 m_adjacency_type;
    int                 m_triplet_type;
    column_vector       m_columns;          //  Graph columns.
    adjacency_vector    m_adjacencies;      //  Graph column adjacencies.
    triplet_vector      m_triplets;         //  Graph column triplets.
    boundary_type       m_boundary;
    point_type          m_origin;

    //  
    // Helper functions
    ///////////////////////////////////////////////////////////////////////
    inline void load_source_xml(const char* filename);
    inline void save_source_xml(std::ostream& os);

    inline void load_target_xml(const char* filename);
    inline void save_target_xml(std::ostream& os);
    
    inline void parse_c(int type, int size, const char* text);
    inline void parse_a(int type, int size, const char* text);
    inline void parse_t(int type, int size, const char* text);
    inline void parse_b(const char* text);
    inline void parse_o(const char* text);
};

} // namespace

#   ifndef __OPTNET_SEPARATION_MODEL__
#       include <optnet/_alpha/graphspec.cxx>
#       if defined(__OPTNET_USE_MSXML__)
#           include <optnet/_alpha/graphspec_msxml.cxx>
#       else
#           include <optnet/_alpha/graphspec_xerces.cxx>        
#       endif
#   endif

#endif // ___GRAPHSPEC_HXX___
