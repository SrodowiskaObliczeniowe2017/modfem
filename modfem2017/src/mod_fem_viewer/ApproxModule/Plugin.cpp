#include "Plugin.h"
#include <iostream>

namespace Approx {

	size_t Plugin::DLLRefCount(0);

	Plugin::hlibModule Plugin::loadPlugin(const std::string & sModuleName)
	{
		hlibModule module_h;
	#ifdef WIN32
		module_h = ::LoadLibrary(sFilename.c_str());
	#else
		// Clear errors
		::dlerror();
		module_h = ::dlopen(sModuleName.c_str(),RTLD_LAZY);
		if(!module_h) {
			std::cerr << "Error!: " << ::dlerror() << std::endl;
		}
	#endif
		return(module_h);
	}

	void Plugin::unloadPlugin(hlibModule pModule)
	{
		if(!pModule) return;

		if(!--DLLRefCount) {
		#ifdef WIN32
			::FreeLibrary(pModule);
		#else
			::dlclose(pModule);
		#endif
		}
	}

	Plugin::Plugin(const std::string & sModuleName) :
     hDLL(0), 
	 pfnGetEngineVersion(0), 
	 pfnRegisterPlugin(0)
	{

		if(!(hDLL = loadPlugin(sModuleName)))
			throw std::runtime_error(std::string("Can't load library: "+sModuleName));

		try {
#ifdef WIN32
			pfnGetEngineVersion = reinterpret_cast<fnGetEngineVersion*>(
				::GetProcAddress(hDLL,"GetEnginVersion"));
			pfnRegisterPlugin = reinterpret_cast<fnRegisterPlugin*>(
				::GetProcAddress(hDLL,"RegisterPlugin"));

			// Check if it is ok
			if(!pfnGetEngineVersion || !pfnRegisterPlugin) {
				throw std::runtime_error(std::string(sModuleName) + " is not valid library!");
			}
#else
			pfnGetEngineVersion = reinterpret_cast<fnGetEngineVersion*>(
							::dlsym(hDLL,"GetEnginVersion"));
			pfnRegisterPlugin = reinterpret_cast<fnRegisterPlugin*>(
							::dlsym(hDLL,"RegisterPlugin"));

			// Check if it is ok
			if(!pfnGetEngineVersion || !pfnRegisterPlugin) {
				throw std::runtime_error(std::string(sModuleName) + " is not valid library! " + ::dlerror());
			}
#endif
			DLLRefCount = 1;
		}
		catch(...)
		{
		#ifdef WIN32
			::FreeLibrary(hDLL);
		#else
			::dlclose(hDLL);
		#endif
			throw(-1);
		}
	}

	Plugin::Plugin(const Plugin& rh) :
		hDLL(rh.hDLL),
		pfnGetEngineVersion(rh.pfnGetEngineVersion),
		pfnRegisterPlugin(rh.pfnRegisterPlugin)
	{
		DLLRefCount++;
	}


	Plugin::~Plugin()
	{
		unloadPlugin(this->hDLL);
	}
}
