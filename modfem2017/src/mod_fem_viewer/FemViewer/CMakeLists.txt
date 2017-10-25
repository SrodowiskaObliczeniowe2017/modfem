# set directories
set(Mod_fem_viewer_dir ${CMAKE_SOURCE_DIR}/mod_fem_viewer)
set(FemViewer_dir ${CMAKE_CURRENT_SOURCE_DIR})
set(FemViewer_incdir "${FemViewer_dir}/inc")
set(FemViewer_srcdir "${FemViewer_dir}/src")
set(FemViewer_resdir "${FemViewer_dir}/res")
set(FemViewer_gltdir "${FemViewer_dir}/glt")
set(FemViewer_windir "${FemViewer_dir}/win")
set(FemViewer_wxdir  "${FemViewer_dir}/wx")
set(FemViewer_ocldir "${FemViewer_dir}/ocl")
set(FemViewer_micdir "${FemViewer_dir}/mic")

# include directory paths
list(APPEND FEM_VIEWER_INCLUDE_DIRS
  ${FemViewer_dir}
  ${FemViewer_incdir}
  ${FemViewer_srcdir}
  ${FemViewer_resdir}
  ${FemViewer_ocldir}
  ${mod_fem_viewer_srcdir}/utils
  ${OPENGL_INCLUDE_DIR}
  ${FREEGLUT_INCLUDE_DIR}
)

# find vtk stuff
if(${MODFEM_MPI} STREQUAL "mpi")
  find_package(VTK REQUIRED)
  include(${VTK_USE_FILE})
endif()
  

# turn on rc complilation
if(WIN32)
  #enable_language(RC)
endif()

# set header files special for VS
unset(HEADER_FILES)
set(HEADER_FILES
  ${CMAKE_SOURCE_DIR}/include/mod_fem_viewer.h
  #${Mod_fem_viewer_dir}/include/mmh_inf.h
  #${Mod_fem_viewer_dir}/include/aph_intf.h
  ${Mod_fem_viewer_dir}/include/ArrayT.hpp
  ${Mod_fem_viewer_dir}/include/BaseHandle.h
  ${Mod_fem_viewer_dir}/include/BaseMesh.h
  ${Mod_fem_viewer_dir}/include/BaseField.h
  ${Mod_fem_viewer_dir}/include/defs.h
  ${Mod_fem_viewer_dir}/include/ElemType.h
  ${Mod_fem_viewer_dir}/include/Enums.h
  ${Mod_fem_viewer_dir}/include/fv_compiler.h
  ${Mod_fem_viewer_dir}/include/fv_config.h
  ${Mod_fem_viewer_dir}/include/Point_T.hh
  ${Mod_fem_viewer_dir}/include/Vector_T.hh
  ${FemViewer_incdir}/WindowFV.h
  ${FemViewer_incdir}/ViewManager.h
  ${FemViewer_incdir}/ModelControler.h
  ${FemViewer_incdir}/Object.h
  ${FemViewer_incdir}/fv_txt_utls.h
  ${FemViewer_incdir}/BBox3D.h
  ${FemViewer_incdir}/GraphicData.h
  ${FemViewer_incdir}/GraphicsSettings.h
  #${FemViewer_incdir}/CubeGraphicsetts
  ${FemViewer_incdir}/Legend.h
  ${FemViewer_incdir}/Tile.h
  ${FemViewer_incdir}/View.h
  ${FemViewer_incdir}/Matrix.h
  ${FemViewer_incdir}/Mesh.h
  ${FemViewer_incdir}/Field.h
  ${FemViewer_incdir}/PrimitiveAccumulator.h
  ${FemViewer_incdir}/Color.h
  #${FemViewer_incdir}/Fire.h
  ${FemViewer_incdir}/Vec3D
  ${FemViewer_incdir}/VtxAccumulator.h
  ${FemViewer_incdir}/VtxPrimitiveAccumulator.h
  ${FemViewer_incdir}/MathHelper.h
  ${FemViewer_incdir}/CutPlane.h
  ${FemViewer_incdir}/Plane.h
  ${FemViewer_incdir}/elem_tables.h
  ${FemViewer_incdir}/Fire.h
  ${FemViewer_incdir}/Log.h
  ${FemViewer_incdir}/ColorBar.h
  ${FemViewer_incdir}/Ray.h
  ${FemViewer_incdir}/Accelerators.h
  ${FemViewer_incdir}/RContext.h
  ${FemViewer_incdir}/Camera.h
  ${FemViewer_incdir}/Geometry.h
  ${FemViewer_incdir}/Shader.h
  ${FemViewer_incdir}/Light.h
  ${FemViewer_incdir}/ExponentEngine.hpp
  ${FemViewer_incdir}/GraphicMesh.hpp
  ${FemViewer_incdir}/GLRenderBase.h
  )


# GUI libraries inlcude dirs etc.
# WINDOWS: WIN32/ATL API
if(WIN32)
  if(GUI_TYPE STREQUAL "WIN")
    list(APPEND FEM_VIEWER_INCLUDE_DIRS ${FemViewer_windir})
    add_definitions(-DVIEWER_WIN_GUI)
    AUX_SOURCE_DIRECTORY(win FV_WIN)
	set(HEADER_FILES ${HEADER_FILES} ${FemViewer_windir}/WinApp.h)
	set(HEADER_FILES ${HEADER_FILES} ${FemViewer_windir}/WinGL.h)
	set(HEADER_FILES ${HEADER_FILES} ${FemViewer_windir}/Handle.h)
	set(HEADER_FILES ${HEADER_FILES} ${FemViewer_windir}/resource.h)
	set(FV_RES ${FemViewer_windir}/mod_fem_viewer_win.rc)
    set_source_files_properties(${FemViewer_windir}/mod_fem_viewer_win.rc PROPERTIES LANGUAGE RC)
	source_group("Resource Files" FILES ${FemViewer_windir}/mod_fem_viewer_win.rc)
    set(FV_GUI_SRC ${FV_WIN})
    message(STATUS "GUI bases on WIN32/ATL library")
    set(VIEWER_GUI_LIBRARY)
  endif()
endif()
 
 # FreeGLUT API
if(GUI_TYPE STREQUAL "GLT")
  list(APPEND FEM_VIEWER_INCLUDE_DIRS ${FemViewer_gltdir})
  add_definitions(-DVIEWER_FGLT_GUI)
  AUX_SOURCE_DIRECTORY(glt FV_GLT)
  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_gltdir}/GlutWindow.h)
  set(FV_GUI_SRC ${FV_GLT})
  message(STATUS "GUI bases on WIN32/ATL library")
  set(VIEWER_GUI_LIBRARY)
 endif()
 
 # search for wxWidgets
if(GUI_TYPE STREQUAL "wxWidgets")
  find_package(wxWidgets COMPONENTS base core gl REQUIRED)
  if(wxWidgets_FOUND)
    if(${GUI_TYPE} STREQUAL "wxWidgets")
      include( ${wxWidgets_USE_FILE})
      if(WIN32)
      	add_definitions(-D__WXMSW_ -D__WIN32__ -DWXUSINGDLL -DwxUSE_GLCANVAS=1 -DwxUSE_GUI=1)
      else(UNIX)
	 add_definitions(-DWXUSINGDLL -DwxUSE_GLCANVAS=1 -DwxUSE_GUI=1)
      endif()
      add_definitions(-DVIEWER_WX_GUI)
      AUX_SOURCE_DIRECTORY(wx  FV_WX)
	  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_windir}/resource.h)
	  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_wxdir}/wxOpenGL.h)
	  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_wxdir}/wxScreenSettings.h)
	  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_wxdir}/wxMainWindow.h)
	  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_wxdir}/wxLegendCanvas.h)
	  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_wxdir}/wxSolCalcDialog.h)
	  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_wxdir}/wxLegendDialog.h)
	  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_wxdir}/wxFemViewerApp.h)
	  set(HEADER_FILES ${HEADER_FILES} ${FemViewer_wxdir}/wxApproxModWnd.h)
	  set(FV_RES ${FemViewer_wxdir}/mod_fem_viewer_wx.rc)
      #set_source_files_properties(${FemViewer_wxdir}/mod_fem_viewer_wx.rc PROPERTIES LANGUAGE RC)
	  source_group("Resource Files" FILES ${FemViewer_wxdir}/mod_fem_viewer_wx.rc)
	  set(FV_GUI_SRC ${FV_WX})
      list(APPEND FEM_VIEWER_INCLUDE_DIRS ${FemViewer_wxdir})
      message(STATUS "GUI bases on wxWidgets libs")
	  message(STATUD " glwx: ${WX_gl_LIBRARY} ")
	  if(WIN32)
        set(VIEWER_GUI_LIBRARY debug wxmsw29d_gl   optimized wxmsw29_gl
                        debug wxmsw29d_core optimized wxmsw29_core
                        debug wxbase29d     optimized wxbase29)
      else(UNIX)
        set(VIEWER_GUI_LIBRARY ${wxWidgets_LIBRARIES})
      endif()
    endif()
  endif()
endif()
 
# set the GUI platform to be used
#set(VIEWER_GUI "wxWidgets")


# set names of target application/library
if(MOD_FEM_VIEWER_USE_AS_LIB)
  set(fvl_dg_prism "fvl_dg_prism")
  set(fvl_dg_prism_2d "fvl_dg_prism2d")
  set(fvl_dg_hybrid "fvl_dg_hybrid")
  set(fvl_std_prism "fvl_std_prism")
  set(fvl_std_prism2d "fvl_std_prism2d")
  set(fvl_std_hybrid "fvl_std_hybrid")
endif()


  
#AUX_SOURCE_DIRECTORY(inc FW_INC)
if(WIN32)
AUX_SOURCE_DIRECTORY(res FW_WIN)
endif()
#file(GLOB_RECURSE INCS "*.h")
source_group("Header Files" FILES ${HEADER_FILES})
unset(SOURCE_FILES)
set(SOURCE_FILES
  ${FemViewer_srcdir}/AppImpl.cpp
  ${FemViewer_srcdir}/mod_fem_viewer
  ${FemViewer_srcdir}/WindowFV
  ${FemViewer_srcdir}/ViewManager
  ${FemViewer_srcdir}/ModelControler
  ${FemViewer_srcdir}/Object
  ${FemViewer_srcdir}/fv_txt_utls
  ${FemViewer_srcdir}/BBox3D
  ${FemViewer_srcdir}/GraphicData
  ${FemViewer_srcdir}/GraphicsSettings
  #${FemViewer_srcdir}/CubeGraphicsetts
  ${FemViewer_srcdir}/Legend
  ${FemViewer_srcdir}/Tile
  ${FemViewer_srcdir}/View
  ${FemViewer_srcdir}/Matrix
  ${FemViewer_srcdir}/Mesh
  ${FemViewer_srcdir}/BaseMesh
  ${FemViewer_srcdir}/Field
  ${FemViewer_srcdir}/BaseField
  ${FemViewer_srcdir}/PrimitiveAccumulator
  ${FemViewer_srcdir}/Color
  #${FemViewer_srcdir}/Fire
  ${FemViewer_srcdir}/Vec3D
  ${FemViewer_srcdir}/VtxAccumulator
  ${FemViewer_srcdir}/VtxPrimitiveAccumulator
  ${FemViewer_srcdir}/MathHelper
  ${FemViewer_srcdir}/CutPlane
  ${FemViewer_srcdir}/Plane
  ${FemViewer_srcdir}/elem_tables
  #${FemViewer_srcdir}/Fire
  ${FemViewer_srcdir}/Log
  ${FemViewer_srcdir}/ColorBar
  ${FemViewer_srcdir}/Ray.cpp
  ${FemViewer_srcdir}/Accelerators.cpp
  ${FemViewer_srcdir}/RContext.cpp
  ${FemViewer_srcdir}/Geometry.cpp
  ${FemViewer_srcdir}/Shader.cpp
  ${FemViewer_srcdir}/Light.cpp
  ${FemViewer_srcdir}/ExponentEngine.cpp
  ${FemViewer_srcdir}/fv_config.cpp
  ${FemViewer_srcdir}/GraphicMesh.cpp
  ${FemViewer_srcdir}/GLRenderBase.cpp
  ${FemViewer_ocldir}/ocl.cpp
  ${FemViewer_micdir}/mic_platform.cpp
  ${FemViewer_micdir}/mic_operator.cpp
  ${FemViewer_micdir}/MicExecutor.cpp
  ${FemViewer_micdir}/MicRenderer.cpp
  ${FV_RES}
)




if(FEM_VIEWER_INCLUDE_DIRS)
  set(FEM_VIEWER_INCLUDE_DIRS ${FEM_VIEWER_INCLUDE_DIRS} ${FEM_VIEWER_EXTERNAL_INCUDE_DIRS} PARENT_SCOPE)
endif()
include_directories(${FEM_VIEWER_INCLUDE_DIRS})

# set additional libraries
set(VIEWER_LIBRARIES fv_utils ${OPENGL_gl_LIBRARY} ${OPENGL_glu_LIBRARY} ${FREEGLUT_LIBRARIES} ${VIEWER_GUI_LIBRARY} 
                      ${GLEW_LIBRARY} 
                      #${VTK_LIBRARIES}
                      )
                      

if(MOD_FEM_VIEWER_USE_AS_LIB)
  add_definitions(-DMOD_FEM_VIEWER)
endif()
    


if(WIN32)
  if(MSVC)
    add_definitions(-D_CONSOLE)		
	# set properties
	#set_target_properties( FemViewer PROPERTIES LINK_FLAGS_DEBUG "/SUBSYSTEM:CONSOLE")
	#set_target_properties( FemViewer PROPERTIES COMPILE_DEFINITIONS_DEBUG "_CONSOLE")
	#set_target_properties( FemViewer PROPERTIES LINK_FLAGS_RELWITHDEBINFO "/SUBSYSTEM:CONSOLE")
	#set_target_properties( FemViewer PROPERTIES COMPILE_DEFINITIONS_RELWITHDEBINFO "_CONSOLE")
	#set_target_properties( FemViewer PROPERTIES LINK_FLAGS_RELEASE "/SUBSYSTEM:CONSOLE")
	#set_target_properties( FemViewer PROPERTIES COMPILE_DEFINITIONS_RELEASE "_CONSOLE")
	#set_target_properties( FemViewer PROPERTIES LINK_FLAGS_MINSIZEREL "/SUBSYSTEM:CONSOLE")
	#set_target_properties( FemViewer PROPERTIES COMPILE_DEFINITIONS_MINSIZEREL "_CONSOLE")
	# create targets	
	if(MOD_FEM_VIEWER_USE_AS_LIB)
	  add_library(${MOD_FEM_VIEWER_LIB} ${SOURCE_FILES} ${FV_GUI_SRC} ${FV_RES} ${HEADER_FILES})
      target_link_libraries(${MOD_FEM_VIEWER_LIB} ${VIEWER_LIBRARIES} C:/Libraries/lib/Win32/Release/libf2c.lib)
	else()
	  add_executable(${MOD_FEM_VIEWER_APP} ${SOURCE_FILES} ${FV_GUI_SRC} ${FW_RES} ${HEADER_FILES})
	  add_dependencies(${MOD_FEM_VIEWER_APP} ApproxModule ApproxDG ApproxSTD ApproxSTDhybrid fv_utils)
   
	
     # if(CMAKE_CFG_INTDIR MATCHES "Debug")
	 #   add_custom_command(TARGET FemViewer POST_BUILD COMMAND "../../../copy32d_dlls.bat")
	 # else(CMAKE_CFG_INTDIR MATCHES "Debug")
	 #   add_custom_command(TARGET FemViewer POST_BUILD COMMAND "../../../copy32_dlls.bat")
	 # endif(CMAKE_CFG_INTDIR MATCHES "Debug")

	

	
      #target_link_libraries(FemViewer IOManager fv_utils ApproxModule ${GLUT_LIB} opengl32 
      #                 debug wxmsw29d_gl   optimized wxmsw29_gl
      #                  debug wxmsw29d_core optimized wxmsw29_core
      #                  debug wxbase29d     optimized wxbase29 )
	  target_link_libraries(${MOD_FEM_VIEWER_LIB} ${VIEWER_LIBRARIES})
	endif()
  endif(MSVC)
elseif(UNIX)
  if(MOD_FEM_VIEWER_USE_AS_LIB)
    add_library(${MOD_FEM_VIEWER_LIB} ${SOURCE_FILES} ${FV_GUI_SRC} ${FV_RES})
    target_link_libraries(${MOD_FEM_VIEWER_LIB} ${VIEWER_LIBRARIES} ${GLUT_LIBRARIES})
  else()   
    add_executable(${MOD_FEM_VIEWER_APP} ${SOURCE_FILES} ${WX_SOURCE_FILES})
    target_link_libraries( ${MOD_FEM_VIEWER_APP} ${MOD_FEM_VIEWER_LIBRARIES} GL GLU ${GLUT_LIBRARIES})
  endif()  
endif(WIN32)

#set(LINK_FLAGS "${LINK_FLAGS} -fopenmp")
SET_TARGET_PROPERTIES(${MOD_FEM_VIEWER_LIB} PROPERTIES COMPILE_DEFINITIONS ${MODFEM_MPI_CODE_DEFINE})
SET_TARGET_PROPERTIES(${MOD_FEM_VIEWER_LIB} PROPERTIES COMPILE_FLAGS "${COMPILE_FLAGS}")
SET_TARGET_PROPERTIES(${MOD_FEM_VIEWER_LIB} PROPERTIES LINK_FLAGS "${LINK_FLAGS}")
#OpenecL kernels
configure_file(${FemViewer_dir}/ocl/kernel_cast_ray.cl ${CMAKE_BINARY_DIR}/kernel_cast_ray.cl COPYONLY)
#configure_file(${FemViewer_dir}/ocl/kernels.cl ${CMAKE_BINARY_DIR}/kernels.cl COPYONLY)
#GLSL shaders
configure_file(${FemViewer_dir}/src/shEdge.frag ${CMAKE_BINARY_DIR}/shEdge.frag COPYONLY)
configure_file(${FemViewer_dir}/src/shEdge.vert ${CMAKE_BINARY_DIR}/shEdge.vert COPYONLY)
configure_file(${FemViewer_dir}/src/shTriVGF.vert ${CMAKE_BINARY_DIR}/shTriVGF.vert COPYONLY)
configure_file(${FemViewer_dir}/src/shTriVGF.frag ${CMAKE_BINARY_DIR}/shTriVGF.frag COPYONLY)
configure_file(${FemViewer_dir}/src/shTriVGF.geom ${CMAKE_BINARY_DIR}/shTriVGF.geom COPYONLY)
configure_file(${FemViewer_dir}/src/shTriStripVGF.geom ${CMAKE_BINARY_DIR}/shTriStripVGF.geom COPYONLY)
