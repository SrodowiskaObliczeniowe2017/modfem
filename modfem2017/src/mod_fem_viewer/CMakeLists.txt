cmake_minimum_required(VERSION 2.8.5)

message("=> Creating graphic module") 
if(MOD_FEM_VIEWER_USE_AS_LIB)
  message(STATUS "=> Graphic module status: external module")
else()
  message(STATUS "=> Graphic module status: standalone application")
endif()

# check if mode_fem_viewer can be build
# search for OpenGL
find_package(OpenGL REQUIRED)
if(OPENGL_FOUND)
  message(STATUS "=> OpenGL found successfully")
  set(BUILD_VIEWER TRUE)
else()
  message(ERROR " OpenGL not found")
  set(BUILD_VIEWER FALSE)
endif()

# search for FreeGLUT
find_package(GLUT REQUIRED)
if(GLUT_FOUND)
  message(STATUS "=> FreeGLUT found successfully")
  set(BUILD_VIEWER TRUE) 
else()
  message(ERROR " FreeGLUT not found!")
  set(BUILD_VIEWER FALSE)
endif()

#VTK for parallel rendering
if(${MODFEM_MPI} STREQUAL "mpi")
  find_package(VTK REQUIRED)
  if(VTK_FOUND)
    message(STATUS "=> VTK found successfully")
    set(BUILD_VIEWER TRUE) 
  else()
    message(ERROR " VTK not found!")
    set(BUILD_VIEWER FALSE)
  endif()
endif()



# search for GLEW


# search for OpenCL
#find_package(OpenCL REQUIRED) [That is doing in main CMAKE file]
#please use opencl configuration in your platform file
#vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
# ------------------------ OPENCL configuration ----------------------- #
# OpenCL - setting opencl acceleration
# Available machine flag: 
# cpu - central processing unit
# gpu - graphics processing unit
# phi - xeon phi acceleartor
#set(OPENCL_INCLUDE_DIRS YOUR_PATH_TO_OPENCL_INCLUDE_DIRECTORY) #Path to directory where is CL/OpenCL folder
#set(OPENCL_LIBRARY_DIRS YOUR_PATH_TO_OPENCL_LIBRARY_DIRECTORY) #Path to directory where is OpenCL library
#^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

include_directories( ${OPENCL_INCLUDE_DIRS} )
if( OPENCL_HAS_CPP_BINDINGS )
message( "OpenCL has CPP bindings. Full include is: " ${OPENCL_INCLUDE_DIRS} )
else( OPENCL_HAS_CPP_BINDINGS )
message( "No OpenCL CPP bindings found" )
endif( OPENCL_HAS_CPP_BINDINGS )
# search components for specifiedtype of GUI
# set defaults:
if(WIN32)
  # set defaults
  set(GUI_TYPE "WIN")
  if(MOD_FEM_VIEWER_GUI_WIN STREQUAL "WX")
    find_package(wxWidgets COMPONENTS base core gl REQUIRED)
	if(wxWidgets_FOUND)
	  message(STATUS "GUI bases on wxWidgets libs")
	  set(GUI_TYPE "wxWidgets")
	else()
	  message(FATAL_ERROR "Can't localize the WxWidghets")
	  set(BUILD_VIEWER FALSE)
	endif()
  elseif(MOD_FEM_VIEWER_GUI_WIN STREQUAL "GLT")
  # seracgh for other type is optional because win32 is allways supported and FreeGLUIT should be found before
    set(GUI_TYPE "GLT")
  endif()
elseif(UNIX)
  # set defaults
  set(GUI_TYPE "wxWidgets")
  if(MOD_FEM_VIEWER_GUI_LIN STREQUAL "WX")
    find_package(wxWidgets)# COMPONENTS base core gl REQUIRED)
    if(wxWidgets_FOUND)
      message(STATUS "GUI bases on wxWidgets libs")
    else()
      message(FATAL_ERROR "Can't localize the WxWidghets")
      set(BUILD_VIEWER FALSE)
    endif()
    
    #Fedora issue
    execute_process(COMMAND wx-config --version OUTPUT_VARIABLE wx_ver)
    if(wx_ver LESS 3.0.0)
      set(WXWIDGETS-3x FALSE)
      message("Fedora issue, set correct WxWidgets version 3.x instead of default 2.x")
      foreach(x RANGE 0 9)
	if(EXISTS ${wxWidgets_CONFIG_EXECUTABLE}-3.${x} AND EXISTS ${wxWidgets_wxrc_EXECUTABLE}-3.${x})

	  #Change executable scripts
	  set(wxWidgets_CONFIG_EXECUTABLE "${wxWidgets_CONFIG_EXECUTABLE}-3.${x}")
	  set(wxWidgets_wxrc_EXECUTABLE "${wxWidgets_wxrc_EXECUTABLE}-3.${x}")

	  #Change include dirs
	  if(EXISTS "/usr/lib64/wx/include/gtk3-unicode-3.${x}")
	    set(wxWidgets_INCLUDE_DIRS "/usr/lib64/wx/include/gtk3-unicode-release-3.${x}")
	  endif()
	  if(EXISTS "/usr/include/wx-3.${x}")
	    set(wxWidgets_INCLUDE_DIRS "${wxWidgets_INCLUDE_DIRS}/usr/include/wx-3.${x}")
	  endif()

	  #Change libraries names
	  foreach(r RANGE 0 9)
	    STRING(REGEX REPLACE "-2.${r}" "-3.${x}" wxWidgets_LIBRARIES ${wxWidgets_LIBRARIES})
	    STRING(REGEX REPLACE "gtk2u" "gtk3u" wxWidgets_LIBRARIES ${wxWidgets_LIBRARIES})
	  endforeach(r)

	  message("New configuration WxWidgets:")
	  message("wxWidgets_CONFIG_EXECUTABLE: ${wxWidgets_CONFIG_EXECUTABLE}")
	  message("wxWidgets_wxrc_EXECUTABLE: ${wxWidgets_wxrc_EXECUTABLE}")
	  message("wxWidgets_INCLUDE_DIRS: ${wxWidgets_INCLUDE_DIRS}")
	  message("wxWidgets_LIBRARIES: ${wxWidgets_LIBRARIES}")
	  set(WXWIDGETS-3x TRUE)
	  break()
	endif()
      endforeach(x)
      if(NOT WXWIDGETS-3x)
	message(FATAL_ERROR "Can't localize the WxWidghets 3.x")
      endif()
    endif()

  elseif(MOD_FEM_VIEWER_GUI_LIN STREQUAL "GLT")
  # seracgh for other type is optional because win32 is allways supported and FreeGLUIT should be found before
    set(GUI_TYPE "GLT")
  endif()
endif()


# TODO:  This involves some erros
#search for GLEW
find_package(GLEW REQUIRED)
if(GLEW_FOUND)
  message("GLEW found successfully")
  #include( ${GLEW_INCLUDE_PATH})
  set(BUILD_VIEWER TRUE)
else()
  set(BUILD_VIEWER FALSE)
endif()

if(BUILD_VIEWER)
  set(my_current_dir ${CMAKE_CURRENT_BINARY_DIR})
  set(my_current_lib_dir "${my_current_dir}/FemViewer")
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${my_current_lib_dir})
  
  set(mod_fem_viewer_srcdir "${CMAKE_SOURCE_DIR}/mod_fem_viewer" CACHE INTERNAL "")
  list(APPEND MOD_FEM_VIEWER_INCLUDE_DIRS
    ${mod_fem_viewer_srcdir}
    ${mod_fem_viewer_srcdir}/include
  )
  if(MOD_FEM_VIEWER_DIRS)
    set(MOD_FEM_VIEWER_DIRS ${MOD_FEM_VIEWER_DIRS} PARENT_SCOPE)
  endif()
  include_directories(${MOD_FEM_VIEWER_INCLUDE_DIRS} ${MOD_FEM_VIEWER_EXTERNAL_INCLUDE_DIRS} ${CMAKE_SOURCE_DIR}/include)
  

  add_subdirectory(utils)
  add_subdirectory(FemViewer)
  if(NOT MOD_FEM_VIEWER_USE_AS_LIB)
    add_subdirectory(ApproxSTDhybrid)
    add_subdirectory(IOManager)
    add_subdirectory(ApproxDG)
    add_subdirectory(ApproxSTD)
    add_subdirectory(ApproxModule) 
  endif()

else()
  message(FATAL_ERROR "=> ERROR: some required components could't be found to build graphic viewer!") 
endif()

