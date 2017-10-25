#ifndef _Accelerators_h__
#define _Accelerators_h__

#include "fv_float.h"
#include "types.h"
#include "Object.h"
#include "Geometry.h"
#include "ElemId.hpp"
#include "MathHelper.h"
#include "BBox3D.h"
#include "Enums.h"
#include "ArrayT.hpp"

#include <vector>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "ocl.h"

namespace FemViewer {

// Forward declarations
//class Object;
class RContext;
class Mesh;
class Field;
struct IsectData;
template<typename TReal,typename TIndex, unsigned N> class mfvObject;
//typedef struct el_isect_into_t;

using namespace fvmath;
template<class t> class Ray;

class Accelerator
{
public:

	Accelerator(const RContext* rc,const char* name_);
	virtual ~Accelerator(void) {}
	virtual const mfvBaseObject* intersect(const Ray<float>& ray, el_isect_info_t *isecData) const;
	virtual void dump(void) const {}
	virtual void drawOn(bool On = true) {}
	virtual void insert(BBox3D&,CVec3f&,ElemId<id_t>) {}
	const RContext *context;
	const char* name;
	Object object;

};

class BVH : public Accelerator
{
	//enum eType { tetra=4, prizm=5, brick=6 };
	static const uint8_t kNumPlaneSetNormals = 3 ;
	static const CVec3f planeSetNormals[kNumPlaneSetNormals];

	struct Extents
	{
		Extents() : object(0), id(-1) {
			for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i)
				d[i][0] = kInfinity, d[i][1] = -kInfinity;
		}
		//~Extents() { elements.clear(); }
		void operator =(const BBox3D& bbox) {
			assert(bbox.isInitialized());
			d[0][0] = bbox.mn.x; d[0][1] = bbox.mx.x;
			d[1][0] = bbox.mn.y; d[1][1] = bbox.mx.y;
			d[2][0] = bbox.mn.z; d[2][1] = bbox.mx.z;
			extbbox = bbox;
		}
		void extendBy(const Extents &extents)
		{
			for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
				if (extents.d[i][0] < d[i][0]) d[i][0] = extents.d[i][0];
				if (extents.d[i][1] > d[i][1]) d[i][1] = extents.d[i][1];
			}
			extbbox += extents.extbbox;
		}
		bool intersect(
			const float *precomputedNumerator,
			const float *precomputeDenominator,
			float &tNear, float &tFar, uint8_t &planeIndex);
		float d[kNumPlaneSetNormals][2]; // d values for each plane-set normals
		const mfvBaseObject *object; // pointer contained by the volume (used by octree)
		int id;
		//std::vector<Elem_Info> elements;
		BBox3D extbbox;
	};
	Extents *extents;
	struct OctreeNode
	{
		OctreeNode *child[8];
		std::vector<const Extents *> data;
		Extents extents;
		bool isLeaf;
		uint8_t depth; // just for debugging
		OctreeNode() : isLeaf(true) {
			memset(child, 0x0, sizeof(OctreeNode *) * 8);
		}
		~OctreeNode() { for (uint8_t i = 0; i < 8; ++i) if (child[i] != NULL) delete child[i]; }
	};
	struct Octree
	{
		Octree(const BBox3D &bbox,Object& obj) : root(NULL),object(obj)
		{
			float sizes[3];
			uint8_t dim = bbox.majorAxis(sizes);
			CVec3f centroid(
				(bbox.mn.x + bbox.mx.x),
				(bbox.mn.y + bbox.mx.y),
				(bbox.mn.z + bbox.mx.z));
			bounds[0] = (centroid - CVec3f(sizes[dim],sizes[dim],sizes[dim])) * 0.5f;
			bounds[1] = (centroid + CVec3f(sizes[dim],sizes[dim],sizes[dim])) * 0.5f;
			root = new OctreeNode;
		}
		void insert(const Extents *extents) { insert(root, extents, bounds[0], bounds[1], 0); }
		void build() { build(root, bounds[0], bounds[1]); }
		void render() { render(root); }
		void drawOn(bool On);
		~Octree() { delete root; }
		struct QueueElement
		{
			const OctreeNode *node; // octree node held by this node in the tree
			float t; // used as key
			QueueElement(const OctreeNode *n, float thit) : node(n), t(thit) {}
    		// comparotor is > instead of < so priority_queue behaves like a min-heap
    		friend bool operator < (const QueueElement &a, const QueueElement &b) { return a.t > b.t; }
		};
		CVec3f bounds[2];
		OctreeNode *root;
		Object& object;
	private:
		void insert(
			OctreeNode *node, const Extents *extents,
			CVec3f boundMin, CVec3f boundMax, int depth)
		{
			if (node->isLeaf) {
				if (node->data.size() == 0 || depth == 16) {
					node->data.push_back(extents);
				}
				else {
					node->isLeaf = false;
					while (node->data.size()) {
						insert(node, node->data.back(), boundMin, boundMax, depth);
						node->data.pop_back();
					}
					insert(node, extents, boundMin, boundMax, depth);
				}
			}
			else {
				// insert bounding volume in the right octree cell
				CVec3f extentsCentroid = (
					CVec3f(extents->d[0][0], extents->d[1][0], extents->d[2][0]) +
					CVec3f(extents->d[0][1], extents->d[1][1], extents->d[2][1])) * 0.5f;
				//mfp_debug("Centroid[%d]: %f %f %f\n",depth,extentsCentroid.x,extentsCentroid.y,extentsCentroid.z);
				CVec3f nodeCentroid = (boundMax + boundMin) * 0.5f;
				CVec3f childBoundMin, childBoundMax;
				uint8_t childIndex = 0;
				if (extentsCentroid[0] > nodeCentroid[0]) {
					childIndex = 4;
					childBoundMin[0] = nodeCentroid[0];
					childBoundMax[0] = boundMax[0];
        		}
        		else {
					childBoundMin[0] = boundMin[0];
					childBoundMax[0] = nodeCentroid[0];
        		}
        		if (extentsCentroid[1] > nodeCentroid[1]) {
					childIndex += 2;
					childBoundMin[1] = nodeCentroid[1];
					childBoundMax[1] = boundMax[1];
        		}
        		else {
					childBoundMin[1] = boundMin[1];
					childBoundMax[1] = nodeCentroid[1];
        		}
        		if (extentsCentroid[2] > nodeCentroid[2]) {
					childIndex += 1;
					childBoundMin[2] = nodeCentroid[2];
					childBoundMax[2] = boundMax[2];
        		}
        		else {
					childBoundMin[2] = boundMin[2];
					childBoundMax[2] = nodeCentroid[2];
        		}
				if (node->child[childIndex] == NULL)
          			node->child[childIndex] = new OctreeNode, node->child[childIndex]->depth = depth;
				insert(node->child[childIndex], extents, childBoundMin, childBoundMax, depth + 1);
			}
		}
		// bottom-up construction
		void build(OctreeNode *node, const CVec3f &boundMin, const CVec3f &boundMax)
		{
			if (node->isLeaf) {
				// compute leaf node bounding volume
				for (uint32_t i = 0; i < node->data.size(); ++i) {
					node->extents.extendBy(*node->data[i]);
				}
			}
			else {
				for (uint8_t i = 0; i < 8; ++i)
					if (node->child[i]) {
						CVec3f pMin, pMax;
						CVec3f pMid = (boundMin + boundMax) * 0.5f;
						pMin[0] = (i & 4) ? pMid[0] : boundMin[0];
						pMax[0] = (i & 4) ? boundMax[0] : pMid[0];
						pMin[1] = (i & 2) ? pMid[1] : boundMin[1];
						pMax[1] = (i & 2) ? boundMax[1] : pMid[1];
						pMin[2] = (i & 1) ? pMid[2] : boundMin[2];
						pMax[2] = (i & 1) ? boundMax[2] : pMid[2];
						build(node->child[i], pMin, pMax);
						node->extents.extendBy(node->child[i]->extents);
					}
			}
		}

		void render(OctreeNode *node/*,uint32_t& totalNodes = 0*/);



	};

	void render() { octree->render(); };
	Octree *octree;
	//const Mesh   *_mesh;
	//const Field  *_field;

public:
	BVH(const RContext *rcx);
	const mfvBaseObject * intersect(const Ray<float>& ray, el_isect_info_t &isectData) const;
	void drawOn(bool On = true) { octree->drawOn(On); }
	~BVH();
};


class Grid : public Accelerator
{
	static const float density;
	static CVec3f calcResolution(const BBox3D& rbox,const uint32_t nelems,uint32_t res[]);

	template<class T>
	struct MyContainer : std::vector<T> {
		CompareBndAct<T> myComp;
		void segregate() {
			std::sort(this->begin(),this->end(), myComp);
		}
	};

	template<class T,class TContainer = MyContainer<T> >
	struct CellItem
	{
		//typedef T value_type;
		//typedef typename T::value_type index_type;

		CellItem() : color(), count(0) {}

		void insert(const ElemInfo &info) {
			vElems.push_back(info);
		}

		void insert(const ElemId<id_t> &ID) {
			vElems.push_back(ID);
		}

		void segregate() {
			vElems.segregate();
		}

		bool intersect(const Ray<float> &ray, const mfvBaseObject **) const;

		TContainer vElems;
		CVec3f color;
		uint32_t count;
	};

	typedef CellItem<ElemId<id_t>, ArrayT<id_t,ElemId<id_t> > > Cell;

public:
	Grid(const RContext *rcx);
	Grid(const RContext* rcx,const BBox3D& bbox,uint32_t nEls);
	~Grid() {
		for (uint32_t i = 0; i < resolution[0] * resolution[1] * resolution[2]; ++i)
			if (cells[i] != NULL) delete cells[i];
		delete [] cells;
		if (C != NULL) delete [] C;
		if (L != NULL) delete [] L;
	}
	Cell* operator()(uint_t x,uint_t y,uint_t z) {
		uint32_t o = z * resolution[0] * resolution[1] + y * resolution[0] + x;
		if (!cells[o]) cells[o] = new Cell;
		return cells[o];
	}
	void test(BBox3D& elBBox,fvmath::CVec3f& center);
	void create();
	void segregate();
	void insert(BBox3D& elBBox,CVec3f& center,ElemId<id_t> ID);
	void PackElems(std::vector<CoordType>& ElCoords,
			       std::vector<float>& ElCoefs,
			       std::vector<uint32_t>& Cells);
	const mfvBaseObject* intersect(const Ray<float> &ray, isect_info_t &isectData) const;
	void render (void);
	void drawOn(bool On = true);
	void dump(void) const;
	Cell **cells;
	BBox3D  bbox;
	const uint32_t  numElems;
	uint32_t resolution[3];
	CVec3f   cellDimension;
	uint32_t* C;
	uint32_t* L;

};

struct gridinfo_t {
	fvmath::Vec3f minDimension;
	fvmath::Vec3f maxDimension;
	fvmath::Vec3f cellSize;
	fvmath::Vec3i cellCount;
};

int create_grid(double targetoccupancy,
		        uint32_t numelems,
		        const BBox3D& gbox,
		        const BBox3D* boxes,
		        cl_int** griddata,
		        cl_int** tridata,
		        gridinfo_t* gridinfo);


}// namespace Femviewer

#endif
