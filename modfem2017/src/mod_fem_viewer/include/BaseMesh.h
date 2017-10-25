#ifndef _BASE_MESH_H_
#define _BASE_MESH_H_

#include <string>

#include "BaseHandle.h"
#include "fv_compiler.h"

namespace FemViewer {


	class BaseMesh : public BaseHandle
	{
	protected:
		/* Pointers to appropriate routines.
		 * When this application is complied and linked as a static module,
		 * routines should be known at once. In other case (standalone program),
		 * it is possible switch between internal and external routines.
		 **/
		static intfint			Get_nmel;
		static intfint			Get_nmno;
		static intf2int			Get_next_act_fa;
		static intf2int			Get_next_act_el;
		static intf2int			Get_face_bc;
		static voidf2int4intp2doublep Get_face_neig;
		static intf2int			Get_el_status;
		static intf2int2intp	Get_el_faces;
		static intf2intintp		Get_el_struct;
		static intf2intintpdoublep Get_el_node_coor;
		static intf2intintp		Get_face_struct;
		static intf2intintpdoublep Get_fa_node_coor;
		static intfint			Get_nmed;
		static intfint			Get_nmfa;
		static intf2intintp 	Get_edge_nodes;
		static intf2int 		Get_edge_status;
		static intf2intdoublep  Get_node_coor;
		static intf2int			Get_node_status;
		static intfcharp		module_introduce;
		static intfintcharpfilep init_mesh;
		static intfint			free_mesh;

		/* Mesh counter
		**/
		static int _meshCounter;
		/* Mesh name or path to input file
		**/
		std::string _name;

	public: // Member routines
		static int GetMeshModuleType();
		/* Constructors
		**/
		explicit BaseMesh(const std::string& name_= "MESH_memory_pool")
		: BaseHandle(), _name(name_) { ; }

		/* Destructor
		**/
		virtual ~BaseMesh() { }

		/* Initialization methods
		**/
		virtual int Init(const char* name_);

		/* Free mesh
		**/
		virtual int Free();

		/* Reload again mesh
		**/
		virtual int Reload(){
			Reset();
			return Init(_name.c_str());
		}

		/* Update internal data
		**/
		virtual int Update(){ return 0;}

		/* Reset to defaults
		 **/
		virtual void Reset() {
			//printf("Basemesh::Reset\n");
			Free();
			BaseHandle::Reset();
		}


		/* Get the Id of the mesh
		 **/
		 	  int& Id() 	  { return this->idx(); }
		const int& Id() const { return this->idx(); }

	    	  std::string& Name() 		{ return _name; }
	    const std::string& Name() const { return _name; }

		/* Specify if mesh is initialized
		 **/
		bool IsInit() const { return ((_meshCounter > 0) && (this->idx() > 0)); }
	};
}


#endif /* _BASE_MESH_H_
*/
