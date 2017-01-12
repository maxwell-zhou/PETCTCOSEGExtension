/*
 ==========================================================================
 |   
 |   $Id: graphspec.cxx 177 2005-02-17 20:26:01Z kangli $
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

#ifndef ___GRAPHSPEC_CXX___
#   define ___GRAPHSPEC_CXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma warning(disable: 4100)
#       pragma warning(disable: 4673)
#       pragma warning(disable: 4671)
#   endif

#   include <optnet/_alpha/graphspec.hxx>
#   include <optnet/_base/except.hxx>
#   include <optnet/_base/secure_s.hxx>
#   include <cstdio>
#   include <set>

namespace optnet {

#   define __GRAPHSPEC_SET_ADJACENCY(a, x, y)   \
        if ((x) < (y)) {                        \
            a.id[0] = (x);                      \
            a.id[1] = (y);                      \
        }                                       \
        else {                                  \
            a.id[0] = (y);                      \
            a.id[1] = (x);                      \
        }

///////////////////////////////////////////////////////////////////////////
void graphspec::adjacencies_from_triplets()
{
    typedef std::set<adjacency_type> set_type;
    
    set_type       edge_set;
    adjacency_type a;

    // Clear current adjacencies.
    m_adjacencies.clear();

    for (triplet_vector::const_iterator it = m_triplets.begin();
         it != m_triplets.end(); ++it)
    {
        __GRAPHSPEC_SET_ADJACENCY(a, it->id[0], it->id[1]);
        if (edge_set.find(a) == edge_set.end()) {
            m_adjacencies.push_back(a);
            edge_set.insert(a);
        }

        __GRAPHSPEC_SET_ADJACENCY(a, it->id[0], it->id[2]);
        if (edge_set.find(a) == edge_set.end()) {
            m_adjacencies.push_back(a);
            edge_set.insert(a);
        }

        __GRAPHSPEC_SET_ADJACENCY(a, it->id[1], it->id[2]);
        if (edge_set.find(a) == edge_set.end()) {
            m_adjacencies.push_back(a);
            edge_set.insert(a);
        }

    }
}

///////////////////////////////////////////////////////////////////////////
void graphspec::positions_from_facts()
{
    double x[OPTNET_GRAPHSPEC_MAX_SURF];
    double y[OPTNET_GRAPHSPEC_MAX_SURF];
    double z[OPTNET_GRAPHSPEC_MAX_SURF];

    for (column_vector::iterator it = m_columns.begin();
         it != m_columns.end(); ++it)
    {
        for (size_type j = 0; j < OPTNET_GRAPHSPEC_MAX_SURF; ++j) {
            x[j] = it->pos[0] + it->normal[0] * it->fact[j];
            y[j] = it->pos[1] + it->normal[1] * it->fact[j];
            z[j] = it->pos[2] + it->normal[2] * it->fact[j];
            it->fact[j] = 0;
        }
        for (size_type j = 0; j < OPTNET_GRAPHSPEC_MAX_SURF; ++j) {
            it->pos[j * 3    ] = x[j];
            it->pos[j * 3 + 1] = y[j];
            it->pos[j * 3 + 2] = z[j];
        }
        //FIXME: A quick hack, should change this later.
        it->normal[3] = it->normal[0];
        it->normal[4] = it->normal[1];
        it->normal[5] = it->normal[2];
    }
}

///////////////////////////////////////////////////////////////////////////
void graphspec::smooth_mesh()
{
    if (m_column_type != 6) {
        throw_exception(std::runtime_error(
            "graphspec::smooth_mesh: Invalid column type."));
    }

    size_type numcols = m_columns.size();
    
    if (numcols == 0) return; // nothing to do.

    size_type numadjs = m_adjacencies.size();

    if (numadjs == 0) {
        // If no adjacencies are present, try to generate them from
        // triplets.
        adjacencies_from_triplets();
    }

    if (numadjs == 0) return; // nothing to do.

    size_type i, j;

    std::vector<std::vector<int> >  column_neighbors(numcols);

    // Create a copy of the columns.
    std::vector<column_type>        columns = m_columns;

    for (i = 0; i < numadjs; ++i) {
        const adjacency_type& adj = m_adjacencies[i];
        column_neighbors[adj.id[0]].push_back(adj.id[1]);
        column_neighbors[adj.id[1]].push_back(adj.id[0]);
    } // i

    for (i = 0; i < numcols; ++i) {
        int n = 1;
        column_type col = columns[i];
        for (j = 0; j < column_neighbors[i].size(); ++j) {
            col += columns[column_neighbors[i][j]];
            ++n;
        } // j
        col /= (double)n;
        m_columns[i] = col;
    } // i
}

///////////////////////////////////////////////////////////////////////////
void graphspec::load_source(const char* filename)
{
    clear(); // Delete all current content.
    m_sorted = false;
    load_source_xml(filename);
    if (!m_sorted) sort_columns();
}

///////////////////////////////////////////////////////////////////////////
void graphspec::save_source(const char* filename)
{
    std::ofstream fs(filename);
    if (fs.is_open()) {
        save_source(fs);
        fs.close();
    }
    else {
        char errmsg[1024];
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::save_source: Could not open file: %s.", filename);
        throw_exception(std::runtime_error(errmsg));
    }
}

///////////////////////////////////////////////////////////////////////////
void graphspec::save_source(std::ostream& os)
{
    if (!m_sorted) sort_columns();
    save_source_xml(os);
}

///////////////////////////////////////////////////////////////////////////
void graphspec::save_source_xml(std::ostream& os)
{
    int sort = m_sorted ? 0 : 1; // need sorting?

    os << "<?xml version=\"1.0\"?>\n";
    os << "<!-- THIS FILE IS AUTOMATICALLY GENERATED. DO NOT EDIT. -->\n";
    os << "<ons version=\"1.0\">\n";
    os << "\t<multicolumn dimension=\"3\">\n";

    // <c> node
    if (m_column_type == 1) {
        if (m_columns.size() > 0) {
            os << "\t\t<c type=\"1\" size=\"" << (int)m_columns.size()
               << "\" sort=\"" << sort << "\">";
            // save all columns
            for (column_vector::const_iterator it1 = m_columns.begin();
                 it1 != m_columns.end();
                 ++it1) {
                os << it1->id     << ','
                   << it1->pos[0] << ','
                   << it1->pos[1] << ','
                   << it1->pos[2] << ';';
            } // for
            os << "</c>\n";
        } else {
            os << "\t\t<c type=\"1\" size=\"0\" sort=\"0\"/>\n";
        }
    }
    else if (m_column_type == 2) {
        if (m_columns.size() > 0) {
            os << "\t\t<c type=\"2\" size=\"" << (int)m_columns.size()
            << "\" sort=\"" << sort << "\">";
            // save all columns
            for (column_vector::const_iterator it1 = m_columns.begin();
                it1 != m_columns.end();
                ++it1) {
                os << it1->id     << ','
                << it1->pos[0]    << ','
                << it1->pos[1]    << ','
                << it1->pos[2]    << ','
                << it1->normal[0] << ','
                << it1->normal[1] << ','
                << it1->normal[2] << ';';
            } // for
            os << "</c>\n";
        } else {
            os << "\t\t<c type=\"2\" size=\"0\" sort=\"0\"/>\n";
        }
    }

    // <a> node
    if (m_adjacencies.size() > 0) {
        os << "\t\t<a type=\"1\" size=\""
           << (int)m_adjacencies.size() << "\">";
        // save all adjacencies
        for (adjacency_vector::const_iterator it2 = m_adjacencies.begin();
             it2 != m_adjacencies.end();
             ++it2) {
            os << it2->id[0] << ','
               << it2->id[1] << ';';
        } // for
        os << "</a>\n";
    } else {
        os << "\t\t<a type=\"1\" size=\"0\"/>\n";
    }
    
    // <t> node
    if (m_triplets.size() > 0) {
        os << "\t\t<t type=\"1\" size=\""
           << (int)m_triplets.size() << "\">";
        // save all triplets
        for (triplet_vector::const_iterator it3 = m_triplets.begin();
             it3 != m_triplets.end();
             ++it3) {
            os << it3->id[0] << ',' 
               << it3->id[1] << ','
               << it3->id[2] << ';';
        } // for
        os << "</t>\n";
    } else {
        os << "\t\t<t type=\"1\" size=\"0\"/>\n";
    }

    os << "\t</multicolumn>\n";
    os << "</ons>\n";
}

///////////////////////////////////////////////////////////////////////////
void graphspec::load_target(const char* filename)
{
    clear(); // Delete all current content.
    m_sorted = false;
    load_target_xml(filename);
    if (!m_sorted) sort_columns();
}

///////////////////////////////////////////////////////////////////////////
void graphspec::save_target(const char* filename)
{
    std::ofstream fs(filename);
    if (fs.is_open()) {
        save_target(fs);
        fs.close();
    }
    else {
        char errmsg[4096];
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::save_target: Could not open file: %s.", filename);
        throw_exception(std::runtime_error(errmsg));
    }
}

///////////////////////////////////////////////////////////////////////////
void graphspec::save_target(std::ostream& os)
{
    if (!m_sorted) sort_columns();
    save_target_xml(os);
}

///////////////////////////////////////////////////////////////////////////
void graphspec::save_target_xml(std::ostream& os)
{
    int sort = m_sorted ? 0 : 1; // need sorting?

    os << "<?xml version=\"1.0\"?>\n";
    os << "<!-- THIS FILE IS AUTOMATICALLY GENERATED. DO NOT EDIT. -->\n";
    os << "<ont version=\"1.0\">\n";
    os << "\t<multicolumn dimension=\"3\" boundary=\""
       << m_boundary[0] << ',' << m_boundary[1] << ',' << m_boundary[2] << ','
       << m_boundary[3] << ',' << m_boundary[4] << ',' << m_boundary[5]
       << "\" origin=\"" 
       << m_origin[0] << ',' << m_origin[1] << ',' << m_origin[2]
       << "\">\n";

    // <c> node
    if (m_column_type == 1) {
        if (m_columns.size() > 0) {
            os << "\t\t<c type=\"3\" size=\"" << (int)m_columns.size()
            << "\" sort=\"" << sort << "\">";
            // save all columns
            for (column_vector::const_iterator it1 = m_columns.begin();
                it1 != m_columns.end();
                ++it1) {
                os << it1->id   << ','
                << it1->pos[0]  << ','
                << it1->pos[1]  << ','
                << it1->pos[2]  << ','
                << it1->fact[0] << ';';
            } // for
            os << "</c>\n";
        } else {
            os << "\t\t<c type=\"3\" size=\"0\" sort=\"0\"/>\n";
        }
    }
    else if (m_column_type == 2 || m_column_type == 4) {
        if (m_columns.size() > 0) {
            os << "\t\t<c type=\"4\" size=\"" << (int)m_columns.size()
            << "\" sort=\"" << sort << "\">";
            // save all columns
            for (column_vector::const_iterator it1 = m_columns.begin();
                it1 != m_columns.end();
                ++it1) {
                os << it1->id     << ','
                << it1->pos[0]    << ','
                << it1->pos[1]    << ','
                << it1->pos[2]    << ','
                << it1->normal[0] << ','
                << it1->normal[1] << ','
                << it1->normal[2] << ','
                << it1->fact[0]   << ';';
            } // for
            os << "</c>\n";
        } else {
            os << "\t\t<c type=\"4\" size=\"0\" sort=\"0\"/>\n";
        }
    }
    else if (m_column_type == 5) {
        if (m_columns.size() > 0) {
            os << "\t\t<c type=\"5\" size=\"" << (int)m_columns.size()
            << "\" sort=\"" << sort << "\">";
            // save all columns
            for (column_vector::const_iterator it1 = m_columns.begin();
                it1 != m_columns.end();
                ++it1) {
                os << it1->id     << ','
                << it1->pos[0]    << ','
                << it1->pos[1]    << ','
                << it1->pos[2]    << ','
                << it1->normal[0] << ','
                << it1->normal[1] << ','
                << it1->normal[2] << ','
                << it1->fact[0]   << ','
                << it1->fact[1]   << ';';
            } // for
            os << "</c>\n";
        } else {
            os << "\t\t<c type=\"5\" size=\"0\" sort=\"0\"/>\n";
        }
    }
    else if (m_column_type == 6) {
        if (m_columns.size() > 0) {
            os << "\t\t<c type=\"6\" size=\"" << (int)m_columns.size()
            << "\" sort=\"" << sort << "\">";
            // save all columns
            for (column_vector::const_iterator it1 = m_columns.begin();
                it1 != m_columns.end();
                ++it1) {
                os << it1->id     << ','
                << it1->pos[0]    << ','
                << it1->pos[1]    << ','
                << it1->pos[2]    << ','
                << it1->normal[0] << ','
                << it1->normal[1] << ','
                << it1->normal[2] << ','
                << it1->pos[3]    << ','
                << it1->pos[4]    << ','
                << it1->pos[5]    << ','
                << it1->normal[3] << ','
                << it1->normal[4] << ','
                << it1->normal[5] << ';';
            } // for
            os << "</c>\n";
        } else {
            os << "\t\t<c type=\"6\" size=\"0\" sort=\"0\"/>\n";
        }
    }
    
    // <t> node
    if (m_triplet_type == 1) {
        if (m_triplets.size() > 0) {
            os << "\t\t<t type=\"1\" size=\""
            << (int)m_triplets.size() << "\">";
            // save all triplets
            for (triplet_vector::const_iterator it2 = m_triplets.begin();
                it2 != m_triplets.end();
                ++it2) {
                os << it2->id[0] << ',' 
                << it2->id[1] << ','
                << it2->id[2] << ';';
            } // for
            os << "</t>\n";
        } else {
            os << "\t\t<t type=\"1\" size=\"0\"/>\n";
        }
    }

    os << "\t</multicolumn>\n";
    os << "</ont>\n";
}

///////////////////////////////////////////////////////////////////////////
void graphspec::parse_c(int type, int size, const char* text)
{
    column_type column;
    char        errmsg[2048];
    char        *token, *context, seps[] = ";";
    int         n = 0;
    
    memset(&column, 0, sizeof(column_type));
    
    if (type == 1) {
        // column format #1: id,x,y,z;
        token = secure_strtok(const_cast<char*>(text), seps, &context);
        while (token != NULL) {
            if (secure_sscanf(token,
                             "%d,%lf,%lf,%lf",
                             &column.id,
                             &(column.pos[0]),
                             &(column.pos[1]),
                             &(column.pos[2])
                ) == 4) {
                
                column.fact[0] = 1.0;
                m_columns.push_back(column);
                if ((0 != size) && (++n >= size)) break;
            }
            else {
                secure_sprintf(errmsg, sizeof(errmsg),
                    "graphspec::parse_c: Failed parsing token: '%s'.",
                    token);
                throw_exception(std::runtime_error(errmsg));
            }
            token = secure_strtok(NULL, seps, &context);
        } // while
    }
    else if (type == 2) {
        // column format #2: id,x,y,z,nx,ny,nz;
        token = secure_strtok(const_cast<char*>(text), seps, &context);
        while (token != NULL) {
            if (secure_sscanf(token,
                              "%d,%lf,%lf,%lf,%lf,%lf,%lf",
                              &column.id,
                              &(column.pos[0]),
                              &(column.pos[1]),
                              &(column.pos[2]),
                              &(column.normal[0]),
                              &(column.normal[1]),
                              &(column.normal[2])
                ) == 7) {
                
                m_columns.push_back(column);
                if ((0 != size) && (++n >= size)) break;
            }
            else {
                secure_sprintf(errmsg, sizeof(errmsg),
                    "graphspec::parse_c: Failed parsing token: '%s'.",
                    token);
                throw_exception(std::runtime_error(errmsg));
            }
            token = secure_strtok(NULL, seps, &context);
        } // while
    }
    else if (type == 3) {
        // column format #3: id,x,y,z,scaling_factor;
        token = secure_strtok(const_cast<char*>(text), seps, &context);
        while (token != NULL) {
            if (secure_sscanf(token,
                             "%d,%lf,%lf,%lf,%lf",
                             &column.id,
                             &(column.pos[0]),
                             &(column.pos[1]),
                             &(column.pos[2]),
                             &(column.fact[0])
                ) == 5) {
                
                m_columns.push_back(column);
                if ((0 != size) && (++n >= size)) break;
            }
            else {
                secure_sprintf(errmsg, sizeof(errmsg),
                    "graphspec::parse_c: Failed parsing token: '%s'.",
                    token);
                throw_exception(std::runtime_error(errmsg));
            }
            token = secure_strtok(NULL, seps, &context);
        } // while
    }
    else if (type == 4) {
        // column format #4: id,x,y,z,nx,ny,nz,scaling_factor;
        token = secure_strtok(const_cast<char*>(text), seps, &context);
        while (token != NULL) {
            if (secure_sscanf(token,
                             "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                             &column.id,
                             &(column.pos[0]),
                             &(column.pos[1]),
                             &(column.pos[2]),
                             &(column.normal[0]),
                             &(column.normal[1]),
                             &(column.normal[2]),
                             &(column.fact[0])
                ) == 8) {
                
                m_columns.push_back(column);
                if ((0 != size) && (++n >= size)) break;
            }
            else {
                secure_sprintf(errmsg, sizeof(errmsg),
                    "graphspec::parse_c: Failed parsing token: '%s'.",
                    token);
                throw_exception(std::runtime_error(errmsg));
            }
            token = secure_strtok(NULL, seps, &context);
        } // while
    }
    else if (type == 5) {
        // column format #5: id,x,y,z,nx,ny,nz,scaling_factor0,scaling_factor1;
        token = secure_strtok(const_cast<char*>(text), seps, &context);
        while (token != NULL) {
            if (secure_sscanf(token,
                             "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                             &column.id,
                             &(column.pos[0]),
                             &(column.pos[1]),
                             &(column.pos[2]),
                             &(column.normal[0]),
                             &(column.normal[1]),
                             &(column.normal[2]),
                             &(column.fact[0]),
                             &(column.fact[1])
                ) == 9) {
                
                m_columns.push_back(column);
                if ((0 != size) && (++n >= size)) break;
            }
            else {
                secure_sprintf(errmsg, sizeof(errmsg),
                    "graphspec::parse_c: Failed parsing token: '%s'.",
                    token);
                throw_exception(std::runtime_error(errmsg));
            }
            token = secure_strtok(NULL, seps, &context);
        } // while
    }
    else if (type == 6) {
        // column format #6: id,x0,y0,z0,nx0,ny0,nz0,x1,y1,z1,nx1,ny1,nz1;
        token = secure_strtok(const_cast<char*>(text), seps, &context);
        while (token != NULL) {
            if (secure_sscanf(token,
                             "%d,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf,%lf",
                             &column.id,
                             &(column.pos[0]),
                             &(column.pos[1]),
                             &(column.pos[2]),
                             &(column.normal[0]),
                             &(column.normal[1]),
                             &(column.normal[2]),
                             &(column.pos[3]),
                             &(column.pos[4]),
                             &(column.pos[5]),
                             &(column.normal[3]),
                             &(column.normal[4]),
                             &(column.normal[5])
                ) == 13) {
                
                m_columns.push_back(column);
                if ((0 != size) && (++n >= size)) break;
            }
            else {
                secure_sprintf(errmsg, sizeof(errmsg),
                    "graphspec::parse_c: Failed parsing token: '%s'.",
                    token);
                throw_exception(std::runtime_error(errmsg));
            }
            token = secure_strtok(NULL, seps, &context);
        } // while
    }
    else {
        throw_exception(std::runtime_error(
            "graphspec::parse_c: Unsupported column type."
            ));
    }

    m_column_type = type;
}

///////////////////////////////////////////////////////////////////////////
void graphspec::parse_a(int type, int size, const char* text)
{
    adjacency_type adjacency;
    char           errmsg[1024];
    char           *token, *context, seps[] = ";";
    int            n = 0;

    if (type != 1) {
        throw_exception(std::runtime_error(
            "graphspec::parse_a: Unsupported adjacency type."
            ));
    }

    token = secure_strtok(const_cast<char*>(text), seps, &context);
    while (token != NULL) {
        if (secure_sscanf(token,
                         "%d,%d",
                         &(adjacency.id[0]),
                         &(adjacency.id[1])
            ) == 2) {

            m_adjacencies.push_back(adjacency); 
            if ((0 != size) && (++n >= size)) break;
        }
        else {
            secure_sprintf(errmsg, sizeof(errmsg),
                "graphspec::parse_a: Failed parsing token: '%s'.",
                token);
            throw_exception(std::runtime_error(errmsg));
        }
        token = secure_strtok(NULL, seps, &context);
    } // while

    m_adjacency_type = type;
}

///////////////////////////////////////////////////////////////////////////
void graphspec::parse_t(int type, int size, const char* text)
{
    triplet_type triplet;
    char         errmsg[1024];
    char         *token, *context, seps[] = ";";
    int          n = 0;

    if (type != 1) {
        throw_exception(std::runtime_error(
            "graphspec::parse_t: Unsupported triplet type."
            ));
    }

    token = secure_strtok(const_cast<char*>(text), seps, &context);
    while (token != NULL) {
        if (secure_sscanf(token,
                          "%d,%d,%d",
                          &(triplet.id[0]),
                          &(triplet.id[1]),
                          &(triplet.id[2])
            ) == 3) {

            m_triplets.push_back(triplet); 
            if ((0 != size) && (++n >= size)) break;
        }
        else {
            secure_sprintf(errmsg, sizeof(errmsg),
                "graphspec::parse_t: Failed parsing token: '%s'.",
                token);
            throw_exception(std::runtime_error(errmsg));
        }
        token = secure_strtok(NULL, seps, &context);
    } // while

    m_triplet_type = type;
}

///////////////////////////////////////////////////////////////////////////
void graphspec::parse_b(const char* text)
{
    if (secure_sscanf(text,
                      "%lf,%lf,%lf,%lf,%lf,%lf",
                      &(m_boundary.v[0]),
                      &(m_boundary.v[1]),
                      &(m_boundary.v[2]),
                      &(m_boundary.v[3]),
                      &(m_boundary.v[4]),
                      &(m_boundary.v[5])) != 6)
    {
        char errmsg[1024];
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::parse_b: Failed parsing token: '%s'.",
            text);
        throw_exception(std::runtime_error(errmsg));
    }
}

///////////////////////////////////////////////////////////////////////////
void graphspec::parse_o(const char* text)
{
    if (secure_sscanf(text,
                      "%lf,%lf,%lf",
                      &(m_origin.v[0]),
                      &(m_origin.v[1]),
                      &(m_origin.v[2])) != 3)
    {
        char errmsg[1024];
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::parse_o: Failed parsing token: '%s'.",
            text);
        throw_exception(std::runtime_error(errmsg));
    }
}

} // namespace

#endif // ___GRAPHSPEC_XERCES_CXX___
