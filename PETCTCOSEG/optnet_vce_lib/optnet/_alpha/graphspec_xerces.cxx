/*
 ==========================================================================
 |   
 |   $Id: graphspec_xerces.cxx 80 2005-01-22 07:49:12Z Administrator $
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

#ifndef ___GRAPHSPEC_XML_CXX___
#   define ___GRAPHSPEC_XML_CXX___

#   if defined(_MSC_VER) && (_MSC_VER > 1000)
#       pragma once
#       pragma warning(disable: 4571)
#       pragma warning(disable: 4100)
#       pragma warning(disable: 4673)
#       pragma warning(disable: 4671)
#       if defined(_DEBUG)
#           pragma comment(lib, "xerces-c_2D.lib")
#           pragma comment(lib, "xerces-depdom_2D.lib")
#       else // !_DEBUG
#           pragma comment(lib, "xerces-c_2.lib")
#           pragma comment(lib, "xerces-depdom_2.lib")
#       endif
#   endif

#   include <optnet/_base/except.hxx>
#   include <optnet/_base/secure_s.hxx>
#   include <optnet/_alpha/graphspec.hxx>
#   include <xercesc/util/OutOfMemoryException.hpp>
#   include <xercesc/util/PlatformUtils.hpp>
#   include <xercesc/util/XMLString.hpp>
#   include <xercesc/dom/DOM.hpp>
#   include <cassert>
#   include <cstdio>

namespace optnet {

    XERCES_CPP_NAMESPACE_USE

    namespace detail {
        // A not-very-efficient XML string to C-string converter
        class strx
        {
            char*  cs;
            strx(const strx&);
            strx& operator=(const strx&);
        public:
            strx(const XMLCh* s) { cs = XMLString::transcode(s); }
            ~strx() { XMLString::release(&cs); }
            operator const char* () const { return cs; }
        };

        // A not-very-efficient C-string to XML string converter
        class xstr
        {
            XMLCh* us;
            xstr(const xstr&);
            xstr& operator=(const xstr&);
        public:
            xstr(const char*  s) { us = XMLString::transcode(s); }
            ~xstr() { XMLString::release(&us); }
            operator const XMLCh*() const { return us; }
        };
    }

///////////////////////////////////////////////////////////////////////////
void graphspec::load_source_xml(const char* filename)
{
    assert(filename != NULL);

    const size_t    MAX_CHARS = 4095;
    char            errmsg[MAX_CHARS + 1];
    bool            error_occured = false;

    try {
        XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& e1) {
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::load_source_xml: %s",
            (const char*)detail::strx(e1.getMessage()));
        throw_exception(std::runtime_error(errmsg));
    }

    XERCES_CPP_NAMESPACE_QUALIFIER DOMImplementation *impl
        = DOMImplementationRegistry::getDOMImplementation(detail::xstr("LS"));
    XERCES_CPP_NAMESPACE_QUALIFIER DOMBuilder *parser
        = ((DOMImplementationLS*)impl)->
            createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
    XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument *doc = 0;
    XERCES_CPP_NAMESPACE_QUALIFIER DOMNodeList *list = 0;
    XERCES_CPP_NAMESPACE_QUALIFIER DOMElement  *root, *spec, *node;

    try {
        doc = parser->parseURI(filename);
    }
    catch (const OutOfMemoryException&) {
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::load_source_xml: Out of memory.");
        error_occured = true;
    }
    catch (const XMLException& e2) {
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::load_source_xml: %s",
            (const char*)detail::strx(e2.getMessage()));
    }
    catch (const DOMException& e3) {
        XMLCh msg[MAX_CHARS + 1];
        if (DOMImplementation::loadDOMExceptionMsg(e3.code, msg, MAX_CHARS))
            secure_sprintf(errmsg, sizeof(errmsg),
                "graphspec::load_source_xml: %s", 
                (const char*)detail::strx(msg));
        else
            secure_sprintf(errmsg, sizeof(errmsg),
                "graphspec::load_source_xml: DOM error code %d",
                e3.code);
        error_occured = true;
    }
    catch (...) {
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::load_source_xml: Failed parsing file %s.",
            filename);
        error_occured = true;
    }

    if (error_occured) {
        parser->release();
        XMLPlatformUtils::Terminate();
        throw_exception(std::runtime_error(errmsg));
    }

    // Begin parsing...
    try {
        root = doc->getDocumentElement();
        
        // ================================================================
        // <ons> node
        if (0 == root ||
            0 != XMLString::compareString
                (root->getTagName(), detail::xstr("ons"))) {
            throw_exception(std::runtime_error(
                "graphspec::load_source_xml: Node 'ons' expected."));
        }

        if (root->hasAttribute(detail::xstr("version")) &&
            0 == XMLString::compareString(root->getAttribute(
                    detail::xstr("version")),
                    detail::xstr("1.0"))) {
            
            // ============================================================
            // <multicolumn> node
            list = root->getElementsByTagName(detail::xstr("multicolumn"));
            if (0 == list || 0 == list->getLength()) {
                throw_exception(std::runtime_error(
                "graphspec::load_source_xml: Node 'multicolumn' expected."
                ));
            }
            
            spec = static_cast<DOMElement*>(list->item(0));
        
            if (spec->hasAttribute(detail::xstr("dimension")) &&
                0 == XMLString::compareString(spec->getAttribute(
                    detail::xstr("dimension")),
                    detail::xstr("3"))) {
            
                const XMLCh* xs;
                int          type, size, sort;

                // ========================================================
                // <c> node
                list = spec->getElementsByTagName(detail::xstr("c"));
                if (0 == list || 0 == list->getLength()) {
                    throw_exception(std::runtime_error(
                        "graphspec::load_source_xml: Node 'c' expected."));
                }

                node = static_cast<DOMElement*>(list->item(0));

                xs = node->getAttribute(detail::xstr("type"));
                if (0 != XMLString::stringLen(xs))
                    type = atoi((const char*)detail::strx(xs));
                else type = 1;

                xs = node->getAttribute(detail::xstr("size"));
                if (0 != XMLString::stringLen(xs))
                    size = atoi((const char*)detail::strx(xs));
                else size = 0;

                xs = node->getAttribute(detail::xstr("sort"));
                if (0 != XMLString::stringLen(xs))
                    sort = atoi((const char*)detail::strx(xs));
                else sort = 0;
                
                m_sorted = (sort == 0); // need sorting?

                // parse columns data
                parse_c(type, size,
                    (const char*)detail::strx(node->getTextContent()));

                // ========================================================
                // <a> node
                list = spec->getElementsByTagName(detail::xstr("a"));
                if (0 == list || 0 == list->getLength()) {
                    throw_exception(std::runtime_error(
                        "graphspec::load_source_xml: Node 'a' expected."));
                }

                node = static_cast<DOMElement*>(list->item(0));

                xs = node->getAttribute(detail::xstr("type"));
                if (0 != XMLString::stringLen(xs))
                    type = atoi((const char*)detail::strx(xs));
                else type = 1;

                xs = node->getAttribute(detail::xstr("size"));
                if (0 != XMLString::stringLen(xs))
                    size = atoi((const char*)detail::strx(xs));
                else size = 0;

                // parse adjacencies data
                parse_a(type, size,
                    (const char*)detail::strx(node->getTextContent()));

                // ========================================================
                // <t> node
                list = spec->getElementsByTagName(detail::xstr("t"));
                if (0 == list || 0 == list->getLength()) {
                    throw_exception(std::runtime_error(
                        "graphspec::load_source_xml: Node 't' expected."));
                }

                node = static_cast<DOMElement*>(list->item(0));

                xs = node->getAttribute(detail::xstr("type"));
                if (0 != XMLString::stringLen(xs))
                    type = atoi((const char*)detail::strx(xs));
                else type = 1;

                xs = node->getAttribute(detail::xstr("size"));
                if (0 != XMLString::stringLen(xs))
                    size = atoi((const char*)detail::strx(xs));
                else size = 0;

                // parse adjacencies data
                parse_t(type, size,
                    (const char*)detail::strx(node->getTextContent()));

            }

        } // if version="1.0"
    }
    catch (...) { 
        parser->release();
        XMLPlatformUtils::Terminate();
        throw;
    }

    parser->release();
    XMLPlatformUtils::Terminate();
}

///////////////////////////////////////////////////////////////////////////
void graphspec::load_target_xml(const char* filename)
{
    assert(filename != NULL);

    const size_t    MAX_CHARS = 4095;
    char            errmsg[MAX_CHARS + 1];
    bool            error_occured = false;

    try {
        XMLPlatformUtils::Initialize();
    }
    catch (const XMLException& e1) {
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::load_target_xml: %s",
            (const char*)detail::strx(e1.getMessage()));
        throw_exception(std::runtime_error(errmsg));
    }

    XERCES_CPP_NAMESPACE_QUALIFIER DOMImplementation *impl
        = DOMImplementationRegistry::getDOMImplementation(detail::xstr("LS"));
    XERCES_CPP_NAMESPACE_QUALIFIER DOMBuilder *parser
        = ((DOMImplementationLS*)impl)->
            createDOMBuilder(DOMImplementationLS::MODE_SYNCHRONOUS, 0);
    XERCES_CPP_NAMESPACE_QUALIFIER DOMDocument *doc = 0;
    XERCES_CPP_NAMESPACE_QUALIFIER DOMNodeList *list = 0;
    XERCES_CPP_NAMESPACE_QUALIFIER DOMElement  *root, *spec, *node;

    try {
        doc = parser->parseURI(filename);
    }
    catch (const OutOfMemoryException&) {
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::load_target_xml: Out of memory.");
        error_occured = true;
    }
    catch (const XMLException& e2) {
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::load_target_xml: %s",
            (const char*)detail::strx(e2.getMessage()));
    }
    catch (const DOMException& e3) {
        XMLCh msg[MAX_CHARS + 1];
        if (DOMImplementation::loadDOMExceptionMsg(e3.code, msg, MAX_CHARS))
            secure_sprintf(errmsg, sizeof(errmsg),
                "graphspec::load_target_xml: %s", 
                (const char*)detail::strx(msg));
        else
            secure_sprintf(errmsg, sizeof(errmsg),
                "graphspec::load_target_xml: DOM error code %d",
                e3.code);
        error_occured = true;
    }
    catch (...) {
        secure_sprintf(errmsg, sizeof(errmsg),
            "graphspec::load_target_xml: Failed parsing file %s.",
            filename);
        error_occured = true;
    }

    if (error_occured) {
        parser->release();
        XMLPlatformUtils::Terminate();
        throw_exception(std::runtime_error(errmsg));
    }

    // Begin parsing...
    try {
        root = doc->getDocumentElement();
        
        // ================================================================
        // <ont> node
        if (0 == root ||
            0 != XMLString::compareString
                (root->getTagName(), detail::xstr("ont"))) {
            throw_exception(std::runtime_error(
                "graphspec::load_target_xml: Node 'ont' expected."));
        }

        if (root->hasAttribute(detail::xstr("version")) &&
            0 == XMLString::compareString(root->getAttribute(
                    detail::xstr("version")),
                    detail::xstr("1.0"))) {
            
            // ============================================================
            // <multicolumn> node
            list = root->getElementsByTagName(detail::xstr("multicolumn"));
            if (0 == list || 0 == list->getLength()) {
                throw_exception(std::runtime_error(
                "graphspec::load_target_xml: Node 'multicolumn' expected."
                ));
            }
            
            spec = static_cast<DOMElement*>(list->item(0));
        
            if (spec->hasAttribute(detail::xstr("dimension")) &&
                0 == XMLString::compareString(spec->getAttribute(
                    detail::xstr("dimension")),
                    detail::xstr("3"))) {
            
                const XMLCh* xs;
                int          type, size, sort;

                // ========================================================
                // boundary attribute
                if (!spec->hasAttribute(detail::xstr("boundary"))) {
                    throw_exception(std::runtime_error(
                        "graphspec::load_target_xml: Attribute 'boundary' expected."
                        ));
                }
                
                parse_b((const char*)detail::strx(
                            spec->getAttribute(detail::xstr("boundary"))
                            )
                        );

                // ========================================================
                // origin attribute
                if (!spec->hasAttribute(detail::xstr("origin"))) {
                    throw_exception(std::runtime_error(
                        "graphspec::load_target_xml: Attribute 'origin' expected."
                        ));
                }
                
                parse_o((const char*)detail::strx(
                            spec->getAttribute(detail::xstr("origin"))
                            )
                        );

                // ========================================================
                // <c> node
                list = spec->getElementsByTagName(detail::xstr("c"));
                if (0 == list || 0 == list->getLength()) {
                    throw_exception(std::runtime_error(
                        "graphspec::load_target_xml: Node 'c' expected."));
                }

                node = static_cast<DOMElement*>(list->item(0));
                
                type = 1;   //
                size = 0;   // default attributes
                sort = 0;   //

                xs = node->getAttribute(detail::xstr("type"));
                if (0 != XMLString::stringLen(xs))
                    type = atoi((const char*)detail::strx(xs));

                xs = node->getAttribute(detail::xstr("size"));
                if (0 != XMLString::stringLen(xs))
                    size = atoi((const char*)detail::strx(xs));

                xs = node->getAttribute(detail::xstr("sort"));
                if (0 != XMLString::stringLen(xs))
                    sort = atoi((const char*)detail::strx(xs));
                
                m_sorted = (sort == 0); // need sorting?

                // parse columns data
                parse_c(type, size,
                    (const char*)detail::strx(node->getTextContent()));

                // ========================================================
                // <t> node
                list = spec->getElementsByTagName(detail::xstr("t"));
                if (0 == list || 0 == list->getLength()) {
                    throw_exception(std::runtime_error(
                        "graphspec::load_target_xml: Node 't' expected."));
                }

                node = static_cast<DOMElement*>(list->item(0));

                type = 1; // default attributes
                size = 0; //

                xs = node->getAttribute(detail::xstr("type"));
                if (0 != XMLString::stringLen(xs))
                    type = atoi((const char*)detail::strx(xs));

                xs = node->getAttribute(detail::xstr("size"));
                if (0 != XMLString::stringLen(xs))
                    size = atoi((const char*)detail::strx(xs));

                // parse adjacencies data
                parse_t(type, size,
                    (const char*)detail::strx(node->getTextContent()));

            }

        } // if version="1.0"
    }
    catch (...) { 
        parser->release();
        XMLPlatformUtils::Terminate();
        throw;
    }

    parser->release();
    XMLPlatformUtils::Terminate();
}

} // namespace

#endif // ___GRAPHSPEC_XML_CXX___
