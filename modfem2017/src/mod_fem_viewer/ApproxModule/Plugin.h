#ifndef _PLUGIN_H_
#define _PLUGIN_H_

#include "../include/ApproxModule.h"
#include "../include/fv_compiler.h"
#include <string>



namespace Approx {
	
class ApproxManager;


class Plugin {
    static size_t      DLLRefCount;        // Number of references to the DLL
public:

    #ifdef WIN32
    typedef HMODULE	hlibModule;
	#else
    typedef void*   hlibModule;
	#endif

public:
	// Function for loading module
    static hlibModule loadPlugin(const std::string & sModuleName);
    // Function for unloading module
    static void unloadPlugin(hlibModule pModule);

public:
	// Initialize and load plugin
    explicit Plugin(const std::string & sFilename);
    // Copy existing plugin instance
    Plugin(const Plugin & rh);
    // Unload a plugin
    ~Plugin();

public:
	// Query the plugin for its expected engine version
    int GetEngineVersion() const {
		return pfnGetEngineVersion();
    }

    /// Register the plugin to a kernel
    void registerPlugin(ApproxManager &AM) {
		pfnRegisterPlugin(AM);
    }
    
  private:
    // Not implemented
    Plugin &operator =(const Plugin &Other);

    // Signature for the version query function
    typedef int fnGetEngineVersion();

    // Signature for the plugin's registration function
    typedef void fnRegisterPlugin(ApproxManager &);



    hlibModule         hDLL;                // DLL handle
    fnGetEngineVersion *pfnGetEngineVersion; // Version query function
    fnRegisterPlugin   *pfnRegisterPlugin;   // Plugin registration function

};

} // end namespace Approx

#endif /* _PLUGIN_H_
 */
