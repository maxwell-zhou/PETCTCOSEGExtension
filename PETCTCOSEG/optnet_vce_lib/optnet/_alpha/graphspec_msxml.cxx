 
 /*
 ==========================================================================
 |   
 |   $Id: graphspec_msxml.cxx 80 2005-01-22 07:49:12Z Administrator $
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
#   endif

#   import <msxml3.dll> // use MSXML 3.0

#   include <optnet/_alpha/graphspec.hxx>
#   include <optnet/_base/except.hxx>
#   include <optnet/_base/secure_s.hxx>
#   include <cassert>
#   include <cstdlib>

namespace optnet {

///////////////////////////////////////////////////////////////////////////
void graphspec::load_source_xml(const char* filename)
{
    assert(filename != NULL);

    HRESULT                         hr;
    MSXML2::IXMLDOMDocumentPtr      docPtr;
    MSXML2::IXMLDOMElementPtr       docElementPtr;
    MSXML2::IXMLDOMNamedNodeMapPtr  attrMapPtr;
    MSXML2::IXMLDOMNodePtr          specPtr;
    MSXML2::IXMLDOMNodePtr          nodePtr;
    MSXML2::IXMLDOMNodePtr          attrPtr;
    _bstr_t                         text;
    int                             type, size, sort;

#   define GRAPHSPEC_LOAD_CLEANUP                           \
        if (NULL != docPtr) docPtr.Release();               \
        if (NULL != docElementPtr) docElementPtr.Release(); \
        if (NULL != attrMapPtr) attrMapPtr.Release();       \
        if (NULL != specPtr) specPtr.Release();             \
        if (NULL != nodePtr) nodePtr.Release();             \
        if (NULL != attrPtr) attrPtr.Release();             \
        ::CoUninitialize()

    // Initialize COM.
    hr = ::CoInitialize(NULL);
    if (FAILED(hr)) {
        throw_exception(std::runtime_error(
            "graphspec::load_source_xml: Failed initializing COM."
            ));
    }

    try {
        
        // Create an instance of the XMLDOM interface.
        hr = docPtr.CreateInstance("Msxml2.DOMDocument");
        if (FAILED(hr)) {
            // Throw std::runtime_error exception.
            throw_exception(std::runtime_error(
                "graphspec::load_source_xml: Failed creating Msxml2.DOMDocument."
                ));
        }

        docPtr->async = VARIANT_FALSE; // default - true,
        
        if (docPtr->load(filename) != VARIANT_TRUE) {
            char errmsg[1024];
            secure_sprintf(errmsg, sizeof(errmsg),
                "graphspec::load_source_xml: Failed loading xml document.\n%s\n%s, Line %d",
                (LPCSTR)docPtr->parseError->Getreason(),
                (LPCSTR)docPtr->parseError->Geturl(),
                (int)docPtr->parseError->Getline()
                );
            // Throw std::runtime_error exception.
            throw_exception(std::runtime_error(errmsg));
        }

        // ================================================================
        // <ons> node
        docElementPtr = docPtr->documentElement;

        if (_bstr_t(docElementPtr->baseName) != _bstr_t("ons")) {
            // Throw std::runtime_error exception.
            throw_exception(std::runtime_error(
                "graphspec::load_source_xml: Node 'ons' expected."
                ));
        }

        attrMapPtr = docElementPtr->attributes;
        if (NULL != attrMapPtr)
            attrPtr = attrMapPtr->getNamedItem(_bstr_t("version"));

        if ((NULL != attrPtr) &&
            (_bstr_t(attrPtr->text) == _bstr_t("1.0"))) { // version 1.0
            
            // ============================================================
            // <multicolumn> node
            specPtr = docElementPtr->firstChild;
            
            if (_bstr_t(specPtr->baseName) != _bstr_t("multicolumn")) {
                // Throw std::runtime_error exception.
                throw_exception(std::runtime_error(
                    "graphspec::load_source_xml: Node 'multicolumn' expected."
                    ));
            }

            attrMapPtr = specPtr->attributes;
            if (NULL != attrMapPtr)
                attrPtr = attrMapPtr->getNamedItem(_bstr_t("dimension"));
            
            if ((NULL != attrPtr) &&
                (_bstr_t(attrPtr->text) == _bstr_t("3"))) {
                
                // ========================================================
                // <c> node
                nodePtr = specPtr->firstChild;
                if ((NULL == nodePtr) || 
                    _bstr_t(nodePtr->baseName) != _bstr_t("c")) {
                    // Throw std::runtime_error exception.
                    throw_exception(std::runtime_error(
                        "graphspec::load_source_xml: Node 'c' expected."
                        ));
                }

                text = nodePtr->text;
                attrMapPtr = nodePtr->attributes;
                if (NULL != attrMapPtr) {
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("type"));
                    type       = (NULL != attrPtr) ? 
                                    atoi(_bstr_t(attrPtr->text)) : 1;
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("size"));
                    size       = (NULL != attrPtr) ? 
                                    atoi(_bstr_t(attrPtr->text)) : 0;
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("sort"));
                    sort       = (NULL != attrPtr) ? 
                                    atoi(_bstr_t(attrPtr->text)) : 0;
                }
                else {
                    type = 1;
                    size = 0;
                    sort = 0;
                }
                
                m_sorted = (sort == 0); // need sorting?

                // parse columns data
                parse_c(type, size, text);

                // ========================================================
                // <a> node
                nodePtr = nodePtr->nextSibling;
                if ((NULL == nodePtr) || 
                    _bstr_t(nodePtr->baseName) != _bstr_t("a")) {
                    // Throw std::runtime_error exception.
                    throw_exception(std::runtime_error(
                        "graphspec::load_source_xml: Node 'a' expected."
                        ));
                }

                text = nodePtr->text;
                attrMapPtr = nodePtr->attributes;
                if (NULL != attrMapPtr) {
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("type"));
                    type       = (NULL != attrPtr) ? 
                                    atoi(_bstr_t(attrPtr->text)) : 1;
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("size"));
                    size       = (NULL != attrPtr) ? 
                                    atoi(_bstr_t(attrPtr->text)) : 0;
                }
                else {
                    type = 1;
                    size = 0;
                }
                
                // parse adjacencies data
                parse_a(type, size, text);
                
                // ========================================================
                // <t> node
                nodePtr = nodePtr->nextSibling;
                if ((NULL == nodePtr) || 
                    _bstr_t(nodePtr->baseName) != _bstr_t("t")) {
                    // Throw std::runtime_error exception.
                    throw_exception(std::runtime_error(
                        "graphspec::load_source_xml: Node 't' expected."
                        ));
                }

                text = nodePtr->text;
                attrMapPtr = nodePtr->attributes;
                if (NULL != attrMapPtr) {
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("type"));
                    type       = (NULL != attrPtr) ?
                                    atoi(_bstr_t(attrPtr->text)) : 1;
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("size"));
                    size       = (NULL != attrPtr) ?
                                    atoi(_bstr_t(attrPtr->text)) : 0;
                }
                else {
                    type = 1;
                    size = 0;
                }
                
                // parse triplets data
                parse_t(type, size, text);

            } // if version="1.0"
        }

        // Release XMLDOM interfaces and finalize COM.
        GRAPHSPEC_LOAD_CLEANUP;
    }
    catch (...) {

        // Release XMLDOM interfaces and finalize COM.
        GRAPHSPEC_LOAD_CLEANUP;
        throw;
    }
}

///////////////////////////////////////////////////////////////////////////
void graphspec::load_target_xml(const char* filename)
{
    assert(filename != NULL);

    HRESULT                         hr;
    MSXML2::IXMLDOMDocumentPtr      docPtr;
    MSXML2::IXMLDOMElementPtr       docElementPtr;
    MSXML2::IXMLDOMNamedNodeMapPtr  attrMapPtr;
    MSXML2::IXMLDOMNodePtr          specPtr;
    MSXML2::IXMLDOMNodePtr          nodePtr;
    MSXML2::IXMLDOMNodePtr          attrPtr;
    _bstr_t                         text;
    int                             type, size, sort;

#   define GRAPHSPEC_LOAD_CLEANUP                           \
        if (NULL != docPtr) docPtr.Release();               \
        if (NULL != docElementPtr) docElementPtr.Release(); \
        if (NULL != attrMapPtr) attrMapPtr.Release();       \
        if (NULL != specPtr) specPtr.Release();             \
        if (NULL != nodePtr) nodePtr.Release();             \
        if (NULL != attrPtr) attrPtr.Release();             \
        ::CoUninitialize()

    // Initialize COM.
    hr = ::CoInitialize(NULL);
    if (FAILED(hr)) {
        throw_exception(std::runtime_error(
            "graphspec::load_target_xml: Failed initializing COM."
            ));
    }

    try {
        
        // Create an instance of the XMLDOM interface.
        hr = docPtr.CreateInstance("Msxml2.DOMDocument");
        if (FAILED(hr)) {
            // Throw std::runtime_error exception.
            throw_exception(std::runtime_error(
                "graphspec::load_target_xml: Failed creating Msxml2.DOMDocument."
                ));
        }

        docPtr->async = VARIANT_FALSE; // default - true,
        
        if (docPtr->load(filename) != VARIANT_TRUE) {
            char errmsg[1024];
            secure_sprintf(errmsg, sizeof(errmsg),
                "graphspec::load_target_xml: Failed loading xml document.\n%s\n%s, Line %d",
                (LPCSTR)docPtr->parseError->Getreason(),
                (LPCSTR)docPtr->parseError->Geturl(),
                (int)docPtr->parseError->Getline()
                );
            // Throw std::runtime_error exception.
            throw_exception(std::runtime_error(errmsg));
        }

        // ================================================================
        // <ont> node
        docElementPtr = docPtr->documentElement;

        if (_bstr_t(docElementPtr->baseName) != _bstr_t("ont")) {
            // Throw std::runtime_error exception.
            throw_exception(std::runtime_error(
                "graphspec::load_target_xml: Node 'ont' expected."
                ));
        }

        attrMapPtr = docElementPtr->attributes;
        if (NULL != attrMapPtr)
            attrPtr = attrMapPtr->getNamedItem(_bstr_t("version"));

        if ((NULL != attrPtr) &&
            (_bstr_t(attrPtr->text) == _bstr_t("1.0"))) { // version 1.0
            
            // ============================================================
            // <multicolumn> node
            specPtr = docElementPtr->firstChild;
            
            if (_bstr_t(specPtr->baseName) != _bstr_t("multicolumn")) {
                // Throw std::runtime_error exception.
                throw_exception(std::runtime_error(
                    "graphspec::load_target_xml: Node 'multicolumn' expected."
                    ));
            }

            attrMapPtr = specPtr->attributes;
            if (NULL != attrMapPtr)
                attrPtr = attrMapPtr->getNamedItem(_bstr_t("dimension"));
            
            if ((NULL != attrPtr) &&
                (_bstr_t(attrPtr->text) == _bstr_t("3"))) {

                // ========================================================
                // boundary attribute
                attrPtr = attrMapPtr->getNamedItem(_bstr_t("boundary"));
                if (NULL == attrPtr) {
                    throw_exception(std::runtime_error(
                        "graphspec::load_target_xml: Attribute 'boundary' expected."
                        ));
                }

                text = attrPtr->text;
                parse_b(text); // parse boundary
                
                // ========================================================
                // origin attribute
                attrPtr = attrMapPtr->getNamedItem(_bstr_t("origin"));
                if (NULL == attrPtr) {
                    throw_exception(std::runtime_error(
                        "graphspec::load_target_xml: Attribute 'origin' expected."
                        ));
                }

                text = attrPtr->text;
                parse_o(text); // parse origin

                // ========================================================
                // <c> node
                nodePtr = specPtr->firstChild;
                if ((NULL == nodePtr) || 
                    _bstr_t(nodePtr->baseName) != _bstr_t("c")) {
                    // Throw std::runtime_error exception.
                    throw_exception(std::runtime_error(
                        "graphspec::load_target_xml: Node 'c' expected."
                        ));
                }

                text = nodePtr->text;
                attrMapPtr = nodePtr->attributes;
                if (NULL != attrMapPtr) {
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("type"));
                    type       = (NULL != attrPtr) ? 
                                    atoi(_bstr_t(attrPtr->text)) : 1;
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("size"));
                    size       = (NULL != attrPtr) ? 
                                    atoi(_bstr_t(attrPtr->text)) : 0;
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("sort"));
                    sort       = (NULL != attrPtr) ? 
                                    atoi(_bstr_t(attrPtr->text)) : 0;
                }
                else {
                    type = 1;
                    size = 0;
                    sort = 0;
                }
                
                m_sorted = (sort == 0); // need sorting?

                // parse columns data
                parse_c(type, size, text);
                
                // ========================================================
                // <t> node
                nodePtr = nodePtr->nextSibling;
                if ((NULL == nodePtr) || 
                    _bstr_t(nodePtr->baseName) != _bstr_t("t")) {
                    // Throw std::runtime_error exception.
                    throw_exception(std::runtime_error(
                        "graphspec::load_target_xml: Node 't' expected."
                        ));
                }

                text = nodePtr->text;
                attrMapPtr = nodePtr->attributes;
                if (NULL != attrMapPtr) {
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("type"));
                    type       = (NULL != attrPtr) ?
                                    atoi(_bstr_t(attrPtr->text)) : 1;
                    attrPtr    = attrMapPtr->getNamedItem(_bstr_t("size"));
                    size       = (NULL != attrPtr) ?
                                    atoi(_bstr_t(attrPtr->text)) : 0;
                }
                else {
                    type = 1;
                    size = 0;
                }
                
                // parse triplets data
                parse_t(type, size, text);

            } // if version="1.0"
        }

        // Release XMLDOM interfaces and finalize COM.
        GRAPHSPEC_LOAD_CLEANUP;
    }
    catch (...) {

        // Release XMLDOM interfaces and finalize COM.
        GRAPHSPEC_LOAD_CLEANUP;
        throw;
    }
}

} // namespace

#endif // ___GRAPHSPEC_XML_CXX___
