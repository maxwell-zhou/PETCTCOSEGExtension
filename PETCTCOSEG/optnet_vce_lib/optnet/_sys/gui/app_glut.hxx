/*
 ==========================================================================
 |   
 |   $Id: app_glut.hxx 21 2005-01-14 15:52:31Z kangli $
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

#ifndef ___APP_GLUT_HXX___
#   define ___APP_GLUT_HXX___

#include <optnet/config.h>
#ifdef __OPTNET_OS_WINNT__
#   include <windows.h>
#endif
#include <gl/gl.h>
#include <gl/glut.h>

#ifdef __OPTNET_CC_MSC__
#   pragma comment(lib,"opengl32.lib")
#   pragma comment(lib,"glu32.lib")
#   pragma comment(lib,"glut32.lib")
#endif

/// @namespace optnet
namespace optnet {
    /// @namespace optnet::system
    namespace system {
        /// @namespace optnet::system::gui
        namespace gui {

///////////////////////////////////////////////////////////////////////////
/// @class app_glut
/// @brief GLUT application skeleton.
///
/// Below is an example of using the app_glut class template.
/// 
/// \code
///
/// #include <optnet/_sys/gui/app_glut.hxx>
/// #include <cstdlib>
/// 
/// using namespace optnet::system::gui;
/// 
/// // Define my own GLUT application class
/// class myapp : 
///     public app_glut<myapp>
/// {
///     typedef app_glut<myapp> _Base;
/// 
/// public:
/// 
///     myapp() :
///         _Base(CALLBACK_RESHAPE | CALLBACK_DISPLAY | CALLBACK_KEY_DOWN)
///     {
///         // Tells the base class (app_glut) to register reshape, display
///         // and keyboard callback functions.
///     }
/// 
///     bool on_initialize()
///     {
///         // app_glut is a minimal application skeleton to the extent that
///         // it (almost) does not do anything. It even does not create the
///         // main window.
/// 
///         // You have to create the main window by yourself as follows.
///         glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
///         glutInitWindowSize(500, 500);
///         glutCreateWindow("A Simple GLUT Application");
/// 
///         glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
///         glEnable(GL_DEPTH_TEST);
/// 
///         //TODO: Add your initialization code here...
/// 
///         // Must return true or the application will terminate
///         // immediately.
///         return true;
///     }
/// 
///     void on_reshape(int width, int height)
///     {
///         glViewport(0, 0, (GLint)width, (GLint)height);
///         glMatrixMode(GL_PROJECTION);
///         glLoadIdentity();
///         glOrtho(-5.0, 5.0, -5.0, 5.0, -5.0, 5.0);
///         glMatrixMode(GL_MODELVIEW); 
///     }
/// 
///     void on_display()
///     {
///         glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
/// 
///         //TODO: Add your display code here...
/// 
///         glutSwapBuffers();
///     }
/// 
///     void on_key_down(unsigned char key, int x, int y)
///     {
///         // Uncomment the following to obtain key modifiers,
///         // e.g., GLUT_ACTIVE_SHIFT, GLUT_ACTIVE_CTRL, or
///         // GLUT_ACTIVE_ATL.
///         // int modifiers = glutGetModifiers();
/// 
///         switch (key) {
///         case 'q':
///         case 'Q':
///             exit(0);
///             break;
/// 
///         //TODO: Add you key handlers here...
/// 
///         default:
///             ;
///         }
///     }
/// }; // End of myapp
/// 
/// 
/// int main(int argc, char* argv[])
/// {
///     return myapp::main(argc, argv);
/// }
///
/// \endcode
///
///////////////////////////////////////////////////////////////////////////
template <class _App>
class app_glut
{
    // Not implemented.
    app_glut(const app_glut&);
    app_glut& operator=(const app_glut&);

public:
    /// @enum feature
    enum feature
    {
        CALLBACK_NONE             = 0x00000000,
        CALLBACK_RESHAPE          = 0x00000001,
        CALLBACK_DISPLAY          = 0x00000002,
        CALLBACK_IDLE             = 0x00000004,
        CALLBACK_KEY_DOWN         = 0x00000008,
        CALLBACK_KEY_UP           = 0x00000010,
        CALLBACK_SPECIAL_KEY_DOWN = 0x00000020,
        CALLBACK_SPECIAL_KEY_UP   = 0x00000040,
        CALLBACK_MOUSE_CLICK      = 0x00000080,
        CALLBACK_MOTION           = 0x00000100,
        CALLBACK_PASSIVE_MOTION   = 0x00000200,
        CALLBACK_ALL              = 0x0000FFFF
    };

    app_glut() throw() :
        m_features(CALLBACK_ALL)
    {
    }
    
    app_glut(int features) throw() :
        m_features(features)
    {
    }

    static int main(int argc, char** argv)
    {
        _App* papp = _App::get_singleton();

        // initialize GLUT
        glutInit(&argc, argv);
                
        if (papp->process_command_line(argc, argv)) {
            if (!papp->on_initialize()) return 2;
        }
        else {
            papp->usage();
            return 1;
        }

        // Install callback funtion table
        if (papp->m_features & CALLBACK_RESHAPE)            glutReshapeFunc         (callback_reshape);
        if (papp->m_features & CALLBACK_DISPLAY)            glutDisplayFunc         (callback_display);
        if (papp->m_features & CALLBACK_IDLE)               glutIdleFunc            (callback_idle);
        if (papp->m_features & CALLBACK_KEY_DOWN)           glutKeyboardFunc        (callback_key_down);
        if (papp->m_features & CALLBACK_KEY_UP)             glutKeyboardUpFunc      (callback_key_up);
        if (papp->m_features & CALLBACK_SPECIAL_KEY_DOWN)   glutSpecialFunc         (callback_special_key_down);
        if (papp->m_features & CALLBACK_SPECIAL_KEY_UP)     glutSpecialUpFunc       (callback_special_key_up);
        if (papp->m_features & CALLBACK_MOUSE_CLICK)        glutMouseFunc           (callback_mouse_click);
        if (papp->m_features & CALLBACK_MOTION)             glutMotionFunc          (callback_motion);
        if (papp->m_features & CALLBACK_PASSIVE_MOTION)     glutPassiveMotionFunc   (callback_passive_motion);

        if (!papp->postprocess_command_line(argc, argv)) return 3;

        // Main loop
        glutMainLoop();
        
        return 0;
    }
        
    bool process_command_line(int argc, char** argv)
    {
        return true;
    }
    
    bool postprocess_command_line(int argc, char** argv)
    {
        return true;
    }

    bool on_initialize()
    {
        // user should implement this function.
        return false;
    }
    
    void usage()
    {
    }
        
    // callback functions
    void on_reshape(int width, int height)
    {
    }
    
    void on_display()
    {
    }
    
    void on_idle()
    {
    }
    
    void on_key_down(unsigned char key, int x, int y)
    {
    }
    
    void on_key_up(unsigned char key, int x, int y)
    {
    }
    
    void on_special_key_down(int key, int x, int y)
    {
    }
    
    void on_special_key_up(int key, int x, int y)
    {
    }
    
    void on_mouse_click(int button, int state, int x, int y)
    {
    }
    
    void on_motion(int x, int y)
    {
    }
    
    void on_passive_motion(int x, int y)
    {
    }
    
    static _App* get_singleton()
    {
        static _App obj;
        return static_cast<_App*>(&obj);
    }


protected:

    int m_features;

    //////////////////////////////////////////////////////////////////////////
    static void callback_reshape(int width, int height)
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp) {
            papp->on_reshape(width, height);
            papp->on_display();
        }
    }
    
    //////////////////////////////////////////////////////////////////////////
    static void callback_display()
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp)
            papp->on_display();
    }
    
    //////////////////////////////////////////////////////////////////////////
    static void callback_idle()
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp)
            papp->on_idle();
    }

    //////////////////////////////////////////////////////////////////////////
    static void callback_key_down(unsigned char key, int x, int y)
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp)
            papp->on_key_down(key, x, y);
    }

    //////////////////////////////////////////////////////////////////////////
    static void callback_key_up(unsigned char key, int x, int y)
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp)
            papp->on_key_up(key, x, y);
    }

    //////////////////////////////////////////////////////////////////////////
    static void callback_special_key_down(int key, int x, int y)
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp)
            papp->on_special_key_down(key, x, y);
    }

    //////////////////////////////////////////////////////////////////////////
    static void callback_special_key_up(int key, int x, int y)
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp)
            papp->on_special_key_up(key, x, y);
    }

    //////////////////////////////////////////////////////////////////////////
    static void callback_mouse_click(int button, int state, int x, int y)
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp)
            papp->on_mouse_click(button, state, x, y);
    }

    //////////////////////////////////////////////////////////////////////////
    static void callback_motion(int x, int y)
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp)
            papp->on_motion(x, y);
    }

    //////////////////////////////////////////////////////////////////////////
    static void callback_passive_motion(int x, int y)
    {
        _App* papp = _App::get_singleton();
        if (NULL != papp)
            papp->on_passive_motion(x, y);
    }
    
};

        } // namespace
    } // namespace
} // namespace

#endif // ___APP_GLUT_HXX___
