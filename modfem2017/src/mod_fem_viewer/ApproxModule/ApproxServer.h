#ifndef _APPROX_SERVER_H_
#define _APPROX_SERVER_H_

#include<memory>
#include<vector>
#include<string>
#include<stdexcept>

namespace Approx {

	class IApproximator;

	class ApproxServer {
	public:
	    // Engine approximation module
		struct ApproxEngin {
			// Constructor
			IApproximator * Self;
			// Destructor
			virtual ~ApproxEngin() {}
			// Get name of the approximation driver
			virtual const std::string & GetName() const = 0;
			// Create a approximation module
			virtual IApproximator* CreateApproximator() = 0;
			// Get handle to approximation module
			virtual IApproximator* GetApproximator() = 0;
		};

		// Destructor
		~ApproxServer() {
			vecApproxModule::iterator it = ApproxModules.begin();
			while(it != ApproxModules.end()) {
				delete *it;
				++it;
			}
		}

		// Add new approximation engine
		void AddApproxModule(std::auto_ptr<ApproxEngin> AE) {
			ApproxModules.push_back(AE.release());
		}

		// Get a number of modules
		inline size_t GetModuleCount() const {
			return ApproxModules.size();
		}

		// Access a module at index
		ApproxEngin & GetModule(const std::string& sModuleName) {
			const size_t size = GetModuleCount();
			for(size_t i(0);i<size;++i)
				if(ApproxModules[i]->GetName() == sModuleName) {
					return *ApproxModules[i];
				}

			throw std::runtime_error(sModuleName);
		}

	private:
		// A vector approximation module
		typedef std::vector<ApproxEngin *> vecApproxModule;
		vecApproxModule ApproxModules;

	};
}

#endif /* _APPROX_SERVER_H_
*/
