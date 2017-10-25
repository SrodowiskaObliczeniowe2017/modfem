#ifndef _FIELD_H_
#define _FIELD_H_

#include "fv_float.h"
#include "fv_config.h"
#include "../../utils/fv_assert.h"
#include "../../include/BaseField.h"
#include"CutPlane.h"
#include "Mesh.h"

#ifdef _USE_FV_EXT_MOD
#include "../../include/ApproxModule.h"
#endif
#include "GraphicElem.hpp"
#include <utility>
#include <vector>
#include <string>



#define FV_DEF_P 101
#define FV_DEF_MIN_P_STD 1
#define FV_DEF_MAX_P_STD 9
#define FV_DEF_MIN_P_DG 101
#define FV_DEF_MAX_P_DG 909


namespace FemViewer {

	class ModelControler;
	class Mesh;
	class Object;
	class RenderParams;
	class CutPlane;
	class Legend;


	//class Field;

	class  Field : public BaseField
	{
		friend class ModelControler;

	public:
		explicit Field(Mesh& pmesh_,const std::string& name_ = "");
		~Field();

		int  Init(const int parent_id_,const char* file_name_);

		void CalculateMinMaxValues(double& min_value, double& max_value,
				std::vector<double>& vec_values,const int ctrl =1) const;

		size_t CountCoeffs(std::vector<int>* degree_ptr = nullptr) const ;
		size_t PackCoeffs(std::vector<mfvFloat_t>& Coeffs,std::vector<int>& Degree);
		static uint32_t FillArrayWithCoeffs(const Field *pField,CoordType Coeffs[]);
		int  Free();
		void Reset();
		void Clear();

		void Render() const;

		int GetMinimalDegree() const { return m_min_deg; }
		int GetMaximalDegree() const { return m_max_deg; }

		double GetMinimalValue() const { return m_min_val; }
		double GetMaximalValue() const { return m_max_val; }

	//private:

		int RenderDG(std::vector<Vertex>& out_vertices,std::vector<unsigned>& out_indices,
				std::vector<int>& out_counts,std::vector<int64_t>& out_offsets);

		int RenderDGCutted(const std::vector<CutPlane>& vplanes, std::vector<Vertex>& out_vertices,std::vector<unsigned>& out_indices
				//,std::vector<int>& out_counts,std::vector<int64_t>& out_offsets
				);
		int RenderSTD(
				std::vector<Vertex>& out_vertices,
				std::vector<unsigned>& out_indices
			);



	public:
		template<typename T = double>
		int GetElementDofs(int nel,T Dofs[]) const {
			return get_element_dofs(this->idx(),nel,this->_solution.currSolution,Dofs);
		}
		int GetElementDegree(int nel) const {
			return get_el_pdeg(this->idx(), nel, NULL);
		}
		int GetElementBaseType(int nel) const {
			return get_base_type(this->idx(), nel);
		}
		int GetNumberOfShapeFunc(int nel,int* degree) const {
			(void)get_el_pdeg(this->idx(),nel,degree);
			return get_el_pdeg_numshap(this->idx(), nel, degree);
		}
		int GetDegreeVector(int nel, int Pdeg[]) const;
		int CalculateSolution(int control, int pdeg, int base,
							  double* ElCoords,
							  double* Ploc,
				              double* Coeffs,
				              double* Sol,
				              double* DSolx = NULL,
				              double* DSoly = NULL,
				              double* DSolz = NULL) const;

		//static void CollectSolutionCoefs()
		#ifdef FV_DEBUG
		void dumpField () const;
		#endif
	private:
		/// Reference to associated mesh
		const Mesh& m_mesh;
		/// Name of the field or path to field file
		std::string m_name;
		// Min/max pdeg
		int m_min_deg;
		int m_max_deg;
		/// Min/max ranges
		mutable double m_min_val;
		mutable double m_max_val;
		/// Vectors of element's degs and values
		std::vector<int> m_elem_degs;
		std::vector<double> m_elem_values;
	};

template<>
inline int Field::GetElementDofs<float>(int nel,float Dofs[]) const
{
	double lDofs[1024];
	int nr_dofs = get_element_dofs(this->idx(),nel,this->_solution.currSolution,lDofs);
	for (int i(0); i<nr_dofs; ++i) Dofs[i] = static_cast<float>(lDofs[i]);
	return nr_dofs;
}

}// end namespace

#endif /* _FIELD_H_
*/
