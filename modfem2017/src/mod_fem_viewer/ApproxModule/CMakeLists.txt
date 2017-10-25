set( APR_SRC 
	Plugin.h
 	Plugin.cpp 
	ApproxServer.h
	../include/ApproxModule.h
	ApproxModule.cpp
	ApproxManager.h
	ApproxManager.cpp )
add_library(ApproxModule STATIC ${APR_SRC})