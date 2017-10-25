#include "ApproxManager.h"

namespace Approx {

ApproxManager * ApproxManager::pInstance = NULL;

ApproxManager::ApproxManager() 
: LoadedModules(), 
  ServerModules()
{}

ApproxManager & ApproxManager::GetInstance()
{
	if(pInstance == NULL)
		pInstance = new ApproxManager();
	return *pInstance;
}

void ApproxManager::LoadApproxModule(const std::string& sModuleName)
{

	if(LoadedModules.find(sModuleName) == LoadedModules.end())
	{
        LoadedModules.insert(ApproxModuleMap::value_type(
          sModuleName,
          Plugin(sModuleName)
        )).first->second.registerPlugin(*this);
    }
}

} // end namespace Approx
