#ifndef _ModelCtrl_H_
#define _ModelCtrl_H_

//#include "Common.h"
#include "Geometry.h" // priomz_info, etc
#include "defs.h"
#include "ocl.h"
//#include "BaseHandle.h"
#include "Mesh.h"
#include "Field.h"
#include "GraphicsSettings.h"
#include "GraphicMesh.hpp"
#include "RContext.h"
//#include "GraphicElem.hpp"
#include "CutPlane.h"
#include "BBox3D.h"
#include "../utils/fv_exception.h"

// stl, system
#include<string>
//#include<stdexcept>
#include<cassert>
#include<memory>
#include<map>
#include<vector>
#include<bitset>



namespace FemViewer
{

	// Forword declarations
	class  Object;
	class  VtxAccumulator;
	class  ViewManager;
	class  BBox3D;
	struct SolutionData;


	class ModelCtrl //: public mfvBaseObject
	{

		friend ModelCtrl & ModelCtrlInst(void);

	public:

		enum eCommnad {
			INIT 	= 0,
			CLEAR	= 1,
			UPDATE  = 2,
			INIT_MESH  = (UPDATE << 1),
			INIT_FIELD = (UPDATE << 2),
			INIT_ACCEL = (UPDATE << 3),
			INIT_GL	   = (UPDATE << 4),
			INIT_LEGEND = (UPDATE << 5),
			INIT_SOL    = (UPDATE << 6),
			INIT_PLANE = (UPDATE << 7),
			INIT_DISPLAY = (UPDATE << 8),
	    };
		// Model items
		enum eItem {
			enumAll = 0,
			enumMesh = 0,
			enumField,
			enumLegend,
			enumAccelerator,
			enumGL,
			enumGLParams,
			enumCutPlane,
			enumSolution,
			enumDrawMethod,
			enumNumItems
		};

		typedef enum {
			MESH_STATE =0,
			FIELD_STATE,
			COLORMAP_STATE,
			CUT_PLANE_STATE,
			ACCELERATOR_STATE,
			GLPARAMS_STATE,
			DISPLAY_STATE,
			ALL_ITEM_STATES,
			NUM_STATES = ALL_ITEM_STATES,
		} notify_status_t;

	public:

		// Dtr
		~ModelCtrl();
		
		bool Do(const int oper = INIT, const char* path_ = NULL);
		void Reset(bool eraseall = false);
		void Clear();
		void Destroy(/*enderParams& params*/);
		void Draw(RenderParams& pParams);
		void SetMeshChange() { m_meshChange = true; }
		void SetFieldChange() { m_fieldChange = true; }
		void SetLegendChange() { m_legendChange = true; }
		void SetCutPlaneChange() { m_cutPlaneChange = true; }
		void SetCutPlaneDisplay() { m_displayChange = true; }
		void SetRenderParamsChange() { m_glParamsChange = true; }

		void NotifyChange(notify_status_t item);
		std::shared_ptr<Legend> LoadColorMap(const char* path);

		int InitData();
		int InitData2();
	public:
		std::vector<prism_info<mfvFloat_t> > elData;
		std::vector<coeffs_info<mfvFloat_t> > coeffsData;
		grid_t gridData;
		int*       C_ptr;
		int*	   L_ptr;
	protected:
		// Ctr
		ModelCtrl();
		//TArrayPtrs<BaseHandle*,10> _objArray;
		GraphicsSettings* m_psettings;
		SolutionData CurrentSolution;
		Mesh    m_mesh;
		Field	m_field;
		RContext* m_pRC;

		bool	m_meshChange;
		bool	m_fieldChange;
		bool	m_legendChange;
		bool	m_solutionChange;
		bool    m_accelChange;
		bool    m_glChange;
		bool    m_glParamsChange;
		bool 	m_cutPlaneChange;
		bool 	m_displayChange;


		std::vector<CutPlane>  m_currPlanes;
		host_info_id _hosts;
		RenderManager<GLCore>* m_pgm;
#ifdef PARALLEL
		int 	  _procId;
		int		  _nrProcesses;
#endif
		std::map<std::string,Legend> m_colorMaps;
		std::bitset<NUM_STATES> m_status;
	public:
		SolutionData& GetCurrentSolution() { return CurrentSolution; }
		Field* GetCurrentField() { return &m_field;  }
        Mesh* GetCurrentMesh() { return &m_mesh; }
		std::vector<CutPlane>& GetCutPlanes() { return m_currPlanes; }
		const std::vector<CutPlane>& GetCutPlanes() const { return m_currPlanes; }
		RContext* RenderingContext()  { return m_pRC; }
		const RContext* RenderingContext() const { return m_pRC; }

		inline int GetMeshModuleType() const { return Mesh::GetMeshModuleType(); }
		inline int GetApproximationType() const { return Field::GetApproximationType(); }

		const BBox3D& Boundary() const { return m_mesh.GetMeshBBox3D(); }

	private:
		// Block use of them
		ModelCtrl(const ModelCtrl&);
		ModelCtrl& operator=(const ModelCtrl&);

		// Internal routines
		bool InitMesh(const char* file_name_ = NULL);
		bool InitField(const char* file_name_= NULL,
				       const HandleType type_= Unknown);
		bool InitVertexAccumulator();
		bool InitAccelStruct();

		void EraseMesh(const int& idx_);
		void EraseField(const int& idx_);

		template<class TEntity>
	    TEntity* GetEntity(const int idx_,const char* name_);

		template<class TMap, class T>
		inline T* GetTPtr(TMap& map_, const int & id_);
		template<class T>
		int RenderColoredLines();

		void BeforeRender();
		void AfterRender();
		bool SetBoundaryElems();
		int  TestConnection() const;
		int  InitHosts(int numHosts);

	};


	extern ModelCtrl& ModelCtrlInst(void); 

	extern int init_data();



} // end namespace FemViewer
#endif /* _ModelCtrl_CONTROLER_HPP_
*/
