#ifndef _BASE_FIELD_H_
#define _BASE_FIELD_H_

#include<string>
#include"fv_compiler.h"
#include"BaseHandle.h"

namespace FemViewer {

	extern const char FieldType[];

	struct SolutionData
	{
		int nrSolutions;
		int nrEquetions;
		int currSolution;
		std::string sFormula;
		SolutionData() : nrSolutions(0), nrEquetions(0), currSolution(0), sFormula("v0")
		{}
		SolutionData(const SolutionData& rhs)
		{
			nrSolutions = rhs.nrSolutions;
			nrEquetions = rhs.nrEquetions;
			currSolution = rhs.currSolution;
			sFormula = rhs.sFormula;
		}
		bool IsInit() const {return (nrSolutions > 0); }
		bool operator=(const SolutionData& rhs)
		{
			bool result(false);
		    result = (rhs.nrSolutions == nrSolutions) ? false : true;
			result = (rhs.nrEquetions == nrEquetions) ? false : true;
			result = (rhs.currSolution == currSolution) ? false : true;
			result = (rhs.sFormula == sFormula) ? false : true;
			if(result)
			{
				nrSolutions = rhs.nrSolutions;
				nrEquetions = rhs.nrEquetions;
				currSolution = rhs.currSolution;
				sFormula = rhs.sFormula;
			}
			return result;
		}
		void Clear()
		{
			nrSolutions = 0;
			nrEquetions = 0;
			currSolution = 0;
			sFormula = "v0";
		}
	};

	class BaseField : protected BaseHandle
	{
	protected:
		/// Internal routines
		static int_global_solutiom  apr_global_solution;
		static double_el_calc 		apr_elem_calc_3D_mod;
		static intfcharp			module_introduce;
		static intfchar5intcharpfnp init_field;
		static intfint 				free_field;
		static intf2int 			get_base_type;
		static intfint 				get_nr_sol;
		static intfint 				get_nreq;
		static intf2intintp 		get_el_pdeg;
		static intf2intintp			get_el_pdeg_numshap;
		static intf3intdoublep 		get_element_dofs;

		/// Solution Data
		SolutionData _solution;
		bool _isSolutionChanged;
	public:
		static int GetApproximationType();
		// Get field index
			  int& Id() 		 { return(this->idx()); }
		const int& Id() const { return(this->idx()); }
		// Get handle to field type
			  HandleType& GetType() 	   { return this->type(); }
		const HandleType& GetType() const { return this->type(); }
		// Gets the number of components in solution
			  int& GetNreq() 	   { return(_solution.nrEquetions); }
		const int& GetNreq() const { return(_solution.nrEquetions); }
		// Gets the current number of solutions
			  int& GetNrSol()       { return(_solution.nrSolutions); }
		const int& GetNrSol() const { return(_solution.nrSolutions); }
		// Gets the index of the current solution
			  int& GetCurrSol() 	  { return(_solution.currSolution); }
		const int& GetCurrSol() const { return(_solution.currSolution); }
		// Gets function formula
			  std::string& GetFromula() 	  { return(_solution.sFormula); }
		const std::string& GetFromula() const { return(_solution.sFormula); }
		void SetSolution(const SolutionData& rhs) 
		{
			_isSolutionChanged = (_solution = rhs); 	
		}
		      SolutionData& GetSolution()       { return _solution; }
		const SolutionData& GetSolution() const { return _solution; }

		bool IsSolutionChanged() const { return _isSolutionChanged; }

		// Default constructor
		explicit BaseField(const HandleType& type_ = Unknown)
		: BaseHandle(type_)
		, _solution(), _isSolutionChanged(false)
		{ }

		// Public destructor
		virtual ~BaseField() {}

		// Initialization
		virtual int Init(const char* name, const int parent_id_);

		// Release field data
		virtual int Free();

		// Reload operation
		//virtual int Reload();

		// Do update after somthing
		//virtual int Update();

		// Clear all
		virtual void Reset();

		// Is geometry ini
		//virtual bool IsGeometryInit() const = 0;

		// Init geomerty
		//virtual void InitGeometry(void * pObject = NULL, char type = '\0') {};

		// Reset geometry
		//virtual void ResetGeometry(void *pObject = NULL, char type = '\0') = 0;

	};




}
#endif /* _BASE_FIELD_H_
*/
