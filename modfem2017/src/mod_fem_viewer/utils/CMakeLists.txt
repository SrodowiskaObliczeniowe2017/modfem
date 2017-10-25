set( FV_UTILS_SRC 
 	fv_assert.h 
	fv_conv_utils.h
	fv_conv_utils.cpp
	fv_dictstr.h
	fv_dictstr.c
	fv_exception.h
	fv_mathparser.h
	fv_mathparser.cpp
	fv_str_utils.h
	fv_str_utils.cpp
	fv_threads.h
	fv_threads.cpp
	fv_timer.h
	fv_timer.cpp
	fv_helper_routines.h
	fv_helper_routines.c
	)
if(NOT MOD_FEM_VIEWER_USE_AS_LIB)
  list(APPEND FV_UTILS_SRC external_utils.c)
endif() 

add_library(fv_utils ${FV_UTILS_SRC})
