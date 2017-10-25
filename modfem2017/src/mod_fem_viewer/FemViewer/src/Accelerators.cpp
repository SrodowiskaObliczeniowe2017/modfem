/*
 * Accelerators.cpp
 *
 *  Created on: 28 sty 2014
 *      Author: Paweł Macioł
 */
#include "Accelerators.h"
#include "RContext.h"
#include "Geometry.h"
#include "Ray.h"
#include "Mesh.h"
#include "Field.h"
#include "BBox3D.h"
#include "ViewManager.h"
#include "VtxAccumulator.h"
#include "Log.h"
#include "defs.h"
#include <queue>
#include <vector>
#include <typeinfo>
#include <fstream>
#include <iostream>

namespace FemViewer {

uint64_t numRayVolumeTests = 0;

Accelerator::Accelerator(const RContext *rcx, const char* name_)
: context(rcx)
, name(name_)
, object(name_)
{
	//mfp_debug("Accelerator: %s",name);
}
//Accelerator::Accelerator(Mesh *rcx) : context(rcx) {}

const mfvBaseObject* Accelerator::intersect(const Ray<float> &ray, el_isect_info_t *isectData) const
{
	float tClosest(ray.tmax);
	mfvBaseObject *hitObject = NULL;
	for (size_t i = 0; i < context->objects.size(); ++i) {
		__sync_fetch_and_add(&numRayVolumeTests, 1);
		Ray<float> r(ray);
		id_t id();
		const BBox3D& lbbox = mfvBaseObject::parentMeshPtr->GetMeshBBox3D();
		if (lbbox.intersect(r)) {
			isect_info_t isectDataCurrent;
			if (context->objects[i]->intersect(ray,&isectDataCurrent)) {
				if (isectDataCurrent.t < tClosest && isectDataCurrent.t > ray.tmin) {
					*(static_cast<isect_info_t*>(isectData)) = isectDataCurrent;
					hitObject = context->objects[i];
					tClosest = isectDataCurrent.t;
				}
			}
		}
	}

	return hitObject;
}


BVH::BVH(const RContext *rcx)
: Accelerator(rcx,"BVH")
, extents(NULL)
, octree(NULL)
//, _mesh(mesh)
//, _field(field)
{
	//mfp_debug("BVH: ctr\n");
	BBox3D bvhbbox,elbbox;
	const uint32_t numEls = context->mesh->GetElems().size();
	extents = new Extents[numEls];

	//mfp_debug("Size of Extens: %u and BBox %u\n",sizeof(Extents),sizeof(BBox3D));
	Mesh::arElConstIter it(context->mesh->GetElems().begin());
	Mesh::arElConstIter it_e(context->mesh->GetElems().end());
	CVec3f vrtx;
	for (uint32_t index(0); it != it_e; ++it,++index)
	{
		// Get elemenent verices
		assert(it->nodes[0] == 6 || it->nodes[0] == 4);
		for (int i = 1; i <= it->nodes[0]; ++i) {
			context->mesh->GetNodeCoor(it->nodes[i],vrtx.v);
			elbbox += vrtx;
		}
		extents[index] = elbbox;
		extents[index].id = index;
		bvhbbox += elbbox;
		elbbox.Reset();
	}

	/*for (uint32_t i = 0; i < mesh->GetElems().size(); ++i) {
		for (uint8_t j = 0; j < kNumPlaneSetNormals; ++j) {
			rcx->objects[i]->computeBounds(planeSetNormals[j], extents[i].d[j][0], extents[i].d[j][1]);
		}
		extents[i].object = rcx->objects[i];

		sceneExtents.extendBy(extents[i]);
	}*/
	// create hierarchy
	octree = new Octree(bvhbbox,this->object);
	for (uint32_t i = 0; i < numEls; ++i) {
		octree->insert(extents + i);
		//std::cout<< "Insert: " << i << " extens\n";
	}
	octree->build();
	//mfp_debug("After build BVH\n");
	octree->render();
}

inline bool BVH::Extents::intersect(
	const float *precomputedNumerator, const float *precomputeDenominator,
	float &tNear, float &tFar, uint8_t &planeIndex)
{
	__sync_fetch_and_add(&numRayVolumeTests, 1);
	for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
		float tn = (d[i][0] - precomputedNumerator[i]) / precomputeDenominator[i];
		float tf = (d[i][1] - precomputedNumerator[i]) / precomputeDenominator[i];
		if (precomputeDenominator[i] < 0) std::swap(tn, tf);
		if (tn > tNear) tNear = tn, planeIndex = i;
		if (tf < tFar) tFar = tf;
		if (tNear > tFar) return false; // test for an early stop
	}

	return true;
}

void BVH::Octree::drawOn(bool On)
{
	object.GetCurrentVtxAccumulator().aRenderEdges = On;
}

const mfvBaseObject* BVH::intersect(const Ray<float>& ray, el_isect_info_t &isectData) const
{
	const mfvBaseObject *hitObject = NULL;
	float precomputedNumerator[BVH::kNumPlaneSetNormals], precomputeDenominator[BVH::kNumPlaneSetNormals];
	for (uint8_t i = 0; i < kNumPlaneSetNormals; ++i) {
		precomputedNumerator[i] = Dot(planeSetNormals[i], ray.orig);
		precomputeDenominator[i] = Dot(planeSetNormals[i], ray.dir);;
	}
#if 0
	float tClosest = ray.tmax;
	for (uint32_t i = 0; i < rc->objects.size(); ++i) {
		__sync_fetch_and_add(&numRayVolumeTests, 1);
		float tNear = -kInfinity, tFar = kInfinity;
		uint8_t planeIndex;
		if (extents[i].intersect(precomputedNumerator, precomputeDenominator, tNear, tFar, planeIndex)) {
			IsectData isectDataCurrent;
			if (rc->objects[i]->intersect(ray, isectDataCurrent)) {
				if (isectDataCurrent.t < tClosest && isectDataCurrent.t > ray.tmin) {
					isectData = isectDataCurrent;
					hitObject = rc->objects[i];
					tClosest = isectDataCurrent.t;
				}
			}
		}
	}
#else
	uint8_t planeIndex = 0;
	float tNear = 0, tFar = ray.tmax;
	if (!octree->root->extents.intersect(precomputedNumerator, precomputeDenominator, tNear, tFar, planeIndex)
		|| tFar < 0 || tNear > ray.tmax)
		return NULL;
	float tMin = tFar;

	std::priority_queue<BVH::Octree::QueueElement> queue;
	queue.push(BVH::Octree::QueueElement(octree->root, 0));
	while(!queue.empty() && queue.top().t < tMin) {
		const OctreeNode *node = queue.top().node;
		queue.pop();
		if (node->isLeaf) {
			el_isect_info_t isectDataCurrent;
			unsigned flag(0);
			for (uint32_t i = 0; i < node->data.size(); ++i) {
				//IsectData isectDataCurrent;
				el_isect_info_t isectDataCurrent;
				uint32_t num = node->data[i]->id;
				const ElemId<id_t> rid(context->mesh->GetElems()[num]);
				bool isIntersection;
				if (rid.is_tetra()) {
					Tetra tetra(rid.eid);
					isIntersection = tetra.intersect(ray,&isectDataCurrent);
				} else {
					Prizm prizm(rid.eid);
					isIntersection = prizm.intersect(ray,&isectDataCurrent);
				}
				if (!isIntersection) continue;
				if (!flag++) {
					isectData = isectDataCurrent;
					continue;
				}
				// the nerarest t
				if (isectData.t > isectDataCurrent.t)
					static_cast<isect_info_t>(isectData) = isectDataCurrent;
				// the farest t
				if (isectData.out.t < isectDataCurrent.out.t)
					isectData.out = isectDataCurrent.out;
				}
/*				if (node->data[i]->object.intersect(ray, isectDataCurrent)) {
					if (isectDataCurrent.t < tMin) {
						tMin = isectDataCurrent.t;
						hitObject = node->data[i]->object;
						isectData = isectDataCurrent;
					}*/

		}
		else {
			for (uint8_t i = 0; i < 8; ++i) {
				if (node->child[i] != NULL) {
					float tNearChild = 0, tFarChild = tFar;
					if (node->child[i]->extents.intersect(precomputedNumerator, precomputeDenominator,
						tNearChild, tFarChild, planeIndex)) {
						float t = (tNearChild < 0 && tFarChild >= 0) ? tFarChild : tNearChild;
						queue.push(BVH::Octree::QueueElement(node->child[i], t));
					}
				}
			}
		}
	}
#endif

	return hitObject;
}


void BVH::Octree::render(OctreeNode* node/*,uint32_t& totalNodes*/)
{
	//mfp_debug("BVH: render\n");
	static CVec3f bboxVertcies[8];
	static uint32_t totalNodes;
	// Init handle to GL object
	//if (! object) object = ViewManagerInst().GetGraphicData().GetRootObject().AddNewObject("Octree");
	// init vertex accumulator
	//mfp_debug("tutaj\n");
	VtxAccumulator& vtxAccum( object.GetCurrentVtxAccumulator());
	//mfp_debug("po\n");
	//vtxAccum.UseVBO() = ViewManagerInst().IsVBOUsed();

	if (node->isLeaf) {
		//mfp_debug("BVH: render leaf %u\n",totalNodes);
		uint8_t num = BBox3D::InitVertices(&(node->extents.extbbox),bboxVertcies);
		for (uint8_t i=0;i<num;++i) object.AddMeshNode(bboxVertcies[i]);
		unsigned int edge[2];
		for (uint8_t i=0;i<12;++i) {
			edge[0] = edges[i].begin() + totalNodes;
			edge[1] = edges[i].end() + totalNodes;
			object.AddEdge(edge);
		}
		totalNodes += num;
		FV_CHECK_ERROR_GL();
	} else {
		for (uint8_t i = 0; i < 8; ++i)
			if (node->child[i]) {
				render(node->child[i]);
				//CVec3f pMin, pMax;
				//CVec3f pMid = (boundMin + boundMax) * 0.5f;
				//pMin[0] = (i & 4) ? pMid[0] : boundMin[0];
				//pMax[0] = (i & 4) ? boundMax[0] : pMid[0];
				//pMin[1] = (i & 2) ? pMid[1] : boundMin[1];
				//pMax[1] = (i & 2) ? boundMax[1] : pMid[1];
				//pMin[2] = (i & 1) ? pMid[2] : boundMin[2];
				//pMax[2] = (i & 1) ? boundMax[2] : pMid[2];
				//build(node->child[i], pMin, pMax);
				//node->extents.extendBy(node->child[i]->extents);
			}

	}

}


BVH::~BVH()
{
	//mfp_debug("BVH: dtr %s\n",this->name);
	if (extents) delete [] extents;
	if (octree != NULL) delete octree;
}

const CVec3f BVH::planeSetNormals[BVH::kNumPlaneSetNormals] = {
	CVec3f(1.f, 0.f, 0.f),
	CVec3f(0.f, 1.f, 0.f),
	CVec3f(0.f, 0.f, 1.f) };
	//CVec3f( sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
	//CVec3f(-sqrtf(3) / 3.f,  sqrtf(3) / 3.f, sqrtf(3) / 3.f),
	//CVec3f(-sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f),
	//CVec3f( sqrtf(3) / 3.f, -sqrtf(3) / 3.f, sqrtf(3) / 3.f) };

#define MAX_X_CELLS 128
#define MAX_NUM_CELLS MAX_X_CELLS * MAX_X_CELLS * MAX_X_CELLS
const float Grid::density = 4.0f;

CVec3f Grid::calcResolution(const BBox3D& rbox,const uint32_t nelems,uint32_t res[])
{
	Vec3f sizes;
#if 0
	uint8_t ndim = rbox.majorAxis(sizes.v);
	const float inv_mxS = 1.0f / sizes.v[ndim];
	const float coeff = 3.0f * powf(nelems, 1.0f / 3.0f);
	const float cellsPreUnitdist = coeff * inv_mxS;
	for (uint8_t i = 0; i < 3; ++i) {
		res[i] = std::floor(sizes.v[i] * cellsPreUnitdist);
		res[i] = std::max(uint32_t(1), std::min(res[i], uint32_t(MAX_X_CELLS)));
	}
#else
	(void)rbox.majorAxis(sizes.v);
	const float inv_V = 1.0f / (sizes.x * sizes.y * sizes.z);
	const float coeff = std::cbrt(density * inv_V * nelems);
	for (uint8_t i = 0; i < 3; ++i) {
		res[i] = std::floor(sizes.v[i] * coeff);
		res[i] = std::max(uint32_t(1), std::min(res[i], uint32_t(MAX_X_CELLS)));
	}
#endif
	assert( res[0] != 0 && res[1] != 0 && res[2] != 0);
	return sizes;
}

Grid::Grid(const RContext *rcx)
: Accelerator(rcx,"Grid")
, cells(NULL)
, bbox(context->mesh->ExtractMeshBBoxFromNodes())
, numElems(context->mesh->GetElems().size())
, cellDimension()
, C(NULL)
, L(NULL)

//Grid::Grid(Mesh *rcx,Field *fp) : Accelerator(rcx), cells(NULL) //, cellMemoryPool(NULL)
{
	// const inverted values for mulitiply
	static const float inv_N[] = { 0.25f, 0.166666667f };
	// compute bound of the scene
	uint32_t totalNumTriangles = 0;
	uint32_t totalNumElements = 0;
	assert(numElems > 0);
	// calculate resolution of the grid
	cellDimension = Grid::calcResolution(this->bbox,this->numElems,this->resolution);
	cellDimension /= CVec3f(resolution[0], resolution[1], resolution[2]);
	//mfp_debug("cellDimensio = %f %f %f\n",cellDimension.x,cellDimension.y,cellDimension.z);
	//mfp_debug("resolutions = %d %d %d\n",resolution[0],resolution[1],resolution[2]);
	// allocate memory
	uint32_t nc = resolution[0] * resolution[1] * resolution[2];
	assert(nc >=1 && nc <= MAX_NUM_CELLS);
	cells = new Grid::Cell* [nc];
	// set all pointers to NULL
	memset(cells, 0x0, sizeof(Grid::Cell*) * nc);

	Vec3d Coords[6];
	int nodes[MMC_MAXELVNO+1];
	//mfp_log_debug("sizeof nodes: %u\n",sizeof(nodes));
	BBox3D cellBBox, elBBox;
	// Insert all the prizm/tetra? in the cells
	int Faces[MMC_MAXELFAC+1],Fa_neig[2],nel(0);
	Mesh::arElemsIter it(context->mesh->GetElems().begin());
	const Mesh::arElConstIter it_e(context->mesh->GetElems().end());
	for (; it != it_e; ++it)
	{
		// Get elemenent verices
		assert(it->nodes[0] == 6 || it->nodes[0] == 4);
		CVec3f vrtx, cent;
		for (int i = 1; i <= it->nodes[0]; ++i) {
			context->mesh->GetNodeCoor(it->nodes[i],vrtx.v);
			cent += vrtx;
			elBBox += vrtx;
		}
		cent *= inv_N[(it->nodes[0] - 4) >> 1];
		Vec3f mind = elBBox[0] - this->bbox[0];
		Vec3f maxd = elBBox[1] - this->bbox[0];
		elBBox[0] = mind / cellDimension;
		elBBox[1] = maxd / cellDimension;
		uint32_t zmin = clamp<uint32_t>(std::floor(elBBox[0].z), 0, resolution[2] - 1);
		uint32_t zmax = clamp<uint32_t>(std::floor(elBBox[1].z), 0, resolution[2] - 1);
		uint32_t ymin = clamp<uint32_t>(std::floor(elBBox[0].y), 0, resolution[1] - 1);
		uint32_t ymax = clamp<uint32_t>(std::floor(elBBox[1].y), 0, resolution[1] - 1);
		uint32_t xmin = clamp<uint32_t>(std::floor(elBBox[0].x), 0, resolution[0] - 1);
		uint32_t xmax = clamp<uint32_t>(std::floor(elBBox[1].x), 0, resolution[0] - 1);
		uint32_t o1 = zmin * resolution[0] * resolution[1] + ymin * resolution[0] + xmin;
		if (((xmax-xmin) == 0) && ((ymax - ymin) == 0) && ((zmax - zmin) == 0)) {
			//uint32_t o1 = zmin * resolution[0] * resolution[1] + ymin * resolution[0] + xmin;
			if (cells[o1] == NULL) cells[o1] = new Cell;
			// Set information that element is whole inside the cell's box
			SET_ACTIVE(it->id);
			cells[o1]->insert(*it);
			//cells[o1]->count++;
			continue;
		}
		CVec3f min(xmin*cellDimension.x,ymin*cellDimension.y,zmin*cellDimension.z);
		cellBBox.mn.z = min.z;
		cellBBox.mx.z = min.z + cellDimension.z;
		// loop over all the cells the triangle overlaps and insert
		for (uint32_t z = zmin; z <= zmax; ++z) {
			cellBBox.mn.y = min.y;
			cellBBox.mx.y = min.y + cellDimension.y;
			for (uint32_t y = ymin; y <= ymax; ++y) {
				cellBBox.mn.x += min.x;
				cellBBox.mx.x += min.x + cellDimension.x;
				for (uint32_t x = xmin; x <= xmax; ++x) {
					if (cellBBox.IsInside(cent)) SET_ACTIVE(it->id);
					uint32_t o = z * resolution[0] * resolution[1] + y * resolution[0] + x;
					if (cells[o] == NULL) cells[o] = new Cell;
					cells[o]->insert(*it);
					//cells[o]->count++;
					cellBBox.mn.x += cellDimension.x;
					cellBBox.mx.x += cellDimension.x;
				}
				cellBBox.mn.y += cellDimension.y;
				cellBBox.mx.y += cellDimension.y;
			}
			cellBBox.mn.z += cellDimension.z;
			cellBBox.mx.z += cellDimension.z;
		}
		elBBox.Reset();
	}
	//mfp_debug("befroe render\n");
	render();
}

Grid::Grid(const RContext* rcx,const BBox3D& bb,uint32_t nEls)
: Accelerator(rcx,"Mesh Grid")
, cells(NULL)
, bbox(bb)
, numElems(nEls)
, C(NULL)
, L(NULL)
{
	//mfp_log_debug("Grid ctr\n");
	assert(numElems > 0);
	// calculate resolution of the grid
	cellDimension = Grid::calcResolution(this->bbox,this->numElems,this->resolution);
	cellDimension /= CVec3f(resolution[0], resolution[1], resolution[2]);
	//mfp_debug("cellDimensio = %f %f %f\n",cellDimension.x,cellDimension.y,cellDimension.z);
	//mfp_debug("resolutions = %d %d %d\n",resolution[0],resolution[1],resolution[2]);
	// allocate memory
	uint32_t nc = resolution[0] * resolution[1] * resolution[2];
	assert(nc >=1 && nc <= MAX_NUM_CELLS);
	cells = new Grid::Cell* [nc];
	// set all pointers to NULL
	memset(cells, 0x0, sizeof(Grid::Cell*) * nc);
	// Zeroed counters
	C = new uint32_t[nc+1];
	memset(C, 0x0, sizeof(uint32_t)*(nc+1));
}

void Grid::test(BBox3D& elBBox,fvmath::CVec3f& center)
{
	//mfp_log_debug("Grid insert\n");
	if (!elBBox.isInitialized()) return;
	static BBox3D cellBBox;
	static CVec3f min;

	Vec3f mind = elBBox[0] - this->bbox.mn;
	Vec3f maxd = elBBox[1] - this->bbox.mx;
	elBBox[0] = mind / cellDimension;
	elBBox[1] = maxd / cellDimension;
	uint32_t zmin = clamp<uint32_t>(std::floor(elBBox[0].z), 0, resolution[2] - 1);
	uint32_t zmax = clamp<uint32_t>(std::floor(elBBox[1].z), 0, resolution[2] - 1);
	uint32_t ymin = clamp<uint32_t>(std::floor(elBBox[0].y), 0, resolution[1] - 1);
	uint32_t ymax = clamp<uint32_t>(std::floor(elBBox[1].y), 0, resolution[1] - 1);
	uint32_t xmin = clamp<uint32_t>(std::floor(elBBox[0].x), 0, resolution[0] - 1);
	uint32_t xmax = clamp<uint32_t>(std::floor(elBBox[1].x), 0, resolution[0] - 1);
	uint32_t o1 = zmin * resolution[0] * resolution[1] + ymin * resolution[0] + xmin;
	if (((xmax-xmin) == 0) && ((ymax - ymin) == 0) && ((zmax - zmin) == 0)) {
		C[o1]++;
		return;
	}
	min.x = xmin*cellDimension.x;
	min.y = ymin*cellDimension.y;
	min.z = zmin*cellDimension.z;
	cellBBox.mn.z = min.z;
	cellBBox.mx.z = min.z + cellDimension.z;
	// loop over all the cells the triangle overlaps and insert
	for (uint32_t z = zmin; z <= zmax; ++z) {
		cellBBox.mn.y = min.y;
		cellBBox.mx.y = min.y + cellDimension.y;
		for (uint32_t y = ymin; y <= ymax; ++y) {
			cellBBox.mn.x += min.x;
			cellBBox.mx.x += min.x + cellDimension.x;
			for (uint32_t x = xmin; x <= xmax; ++x) {
				uint32_t o = z * resolution[0] * resolution[1] + y * resolution[0] + x;
				C[o]++;
				cellBBox.mn.x += cellDimension.x;
				cellBBox.mx.x += cellDimension.x;
			}
			cellBBox.mn.y += cellDimension.y;
			cellBBox.mx.y += cellDimension.y;
		}
		cellBBox.mn.z += cellDimension.z;
		cellBBox.mx.z += cellDimension.z;
	}
}

void Grid::create()
{
	uint32_t grid_size = resolution[0] *resolution[1]* resolution[2];
	for (uint32_t i = 1; i <= grid_size; i++) C[i] += C[i-1];
	assert(C[grid_size] >= numElems);
	printf("C[%u] = %u\n",grid_size,C[grid_size]);
	L = new uint32_t [C[grid_size]];
}

void Grid::segregate()
{

}

void Grid::insert(BBox3D& elBBox,fvmath::CVec3f& center,ElemId<id_t> ID)
{
	//mfp_log_debug("Grid insert\n");
	if (!elBBox.isInitialized()) return;
	static BBox3D cellBBox;
	static CVec3f min;

	Vec3f mind = elBBox[0] - this->bbox.mn;
	Vec3f maxd = elBBox[1] - this->bbox.mx;
	elBBox[0] = mind / cellDimension;
	elBBox[1] = maxd / cellDimension;
	uint32_t zmin = clamp<uint32_t>(std::floor(elBBox[0].z), 0, resolution[2] - 1);
	uint32_t zmax = clamp<uint32_t>(std::floor(elBBox[1].z), 0, resolution[2] - 1);
	uint32_t ymin = clamp<uint32_t>(std::floor(elBBox[0].y), 0, resolution[1] - 1);
	uint32_t ymax = clamp<uint32_t>(std::floor(elBBox[1].y), 0, resolution[1] - 1);
	uint32_t xmin = clamp<uint32_t>(std::floor(elBBox[0].x), 0, resolution[0] - 1);
	uint32_t xmax = clamp<uint32_t>(std::floor(elBBox[1].x), 0, resolution[0] - 1);
	uint32_t o1 = zmin * resolution[0] * resolution[1] + ymin * resolution[0] + xmin;
	/*if (((xmax-xmin) == 0) && ((ymax - ymin) == 0) && ((zmax - zmin) == 0)) {
		//uint32_t o1 = zmin * resolution[0] * resolution[1] + ymin * resolution[0] + xmin;
		if (cells[o1] == NULL) cells[o1] = new Cell;
		// Set information that element is whole inside the cell's box
		//SET_ACTIVE(ID.id);
		cells[o1]->insert(ID);
		C[o1]--;
		return;
	}*/
	//SET_BOUND(ID.id);
	min.x = xmin*cellDimension.x;
	min.y = ymin*cellDimension.y;
	min.z = zmin*cellDimension.z;
	cellBBox.mn.z = min.z;
	cellBBox.mx.z = min.z + cellDimension.z;
	// loop over all the cells the triangle overlaps and insert
	for (uint32_t z = zmin; z <= zmax; ++z) {
		cellBBox.mn.y = min.y;
		cellBBox.mx.y = min.y + cellDimension.y;
		for (uint32_t y = ymin; y <= ymax; ++y) {
			cellBBox.mn.x += min.x;
			cellBBox.mx.x += min.x + cellDimension.x;
			for (uint32_t x = xmin; x <= xmax; ++x) {
				//if (cellBBox.IsInside(center)) SET_ACTIVE(ID.id);
				uint32_t o = z * resolution[0] * resolution[1] + y * resolution[0] + x;
				if (cells[o] == NULL) cells[o] = new Cell;
				cells[o]->insert(ID);
				L[ --C[o] ] = ID.id;
				cellBBox.mn.x += cellDimension.x;
				cellBBox.mx.x += cellDimension.x;
			}
			cellBBox.mn.y += cellDimension.y;
			cellBBox.mx.y += cellDimension.y;
		}
		cellBBox.mn.z += cellDimension.z;
		cellBBox.mx.z += cellDimension.z;
	}
}


template<class T,class TContainer>
bool Grid::CellItem<T,TContainer>::intersect(const Ray<float> &ray, const mfvBaseObject **hitObject) const
{
	el_isect_info_t isecData;
	unsigned flag(0);
	for (size_t i = 0; i < vElems.size(); ++i) {
		Prizm prizm(vElems[i].eid);
		el_isect_info_t currentSecData;
		if (prizm.intersect(ray,&currentSecData)) {
			if (!flag++) {
				isecData = currentSecData;
				continue;
			}
			// the nerarest t
			if (isecData.t > currentSecData.t)
				static_cast<isect_info_t>(isecData)  = currentSecData;
			// the farest t
			if (isecData.out.t < currentSecData.out.t)
				isecData.out = currentSecData.out;
		}

	}

	return (flag!=0);
}

void Grid::PackElems(std::vector<CoordType>& ElCoords,
					 std::vector<float>& ElCoefs,
		             std::vector<uint32_t>& Cells)
{
	// Loop over cells
	const uint32_t nc = resolution[0]*resolution[1]*resolution[2];
	for (uint32_t index = 0; index < nc; ++index) {

	}
	/*for (uint32_t z
				for (uint32_t y = ymin; y <= ymax; ++y) {
					for (uint32_t x = xmin; x <= xmax; ++x) {
						uint32_t o = z * resolution[0] * resolution[1] + y * resolution[0] + x;
						if (cells[o] == NULL) cells[o] = new Cell;
						cells[o]->insert(*it);
					}
				}
			}*/
}

const mfvBaseObject* Grid::intersect(const Ray<float> &ray, isect_info_t &isectData) const
{
	// if the ray doesn't intersect the grid return
	Ray<float> r(ray);
	if (!bbox.intersect(r)) return NULL;
	// initialization step
	CVec3i exit, step, cell;
	CVec3f deltaT, nextCrossingT;
	for (uint8_t i = 0; i < 3; ++i) {
		// convert ray starting point to cell coordinates
		float rayOrigCell = ((r.orig[i] + r.dir[i] * r.tmin) -  bbox[0].v[i]);
		cell[i] = clamp<uint32_t>(std::floor(rayOrigCell / cellDimension[i]), 0, resolution[i] - 1);
		if (r.dir[i] < 0) {
			deltaT[i] = -cellDimension[i] * r.invdir[i];
			nextCrossingT[i] = r.tmin + (cell[i] * cellDimension[i] - rayOrigCell) * r.invdir[i];
			exit[i] = -1;
			step[i] = -1;
		}
		else {
			deltaT[i] = cellDimension[i] * r.invdir[i];
			nextCrossingT[i] = r.tmin + ((cell[i] + 1)  * cellDimension[i] - rayOrigCell) * r.invdir[i];
			exit[i] = resolution[i];
			step[i] = 1;
		}
	}

	// walk through each cell of the grid and test for an intersection if
	// current cell contains geometry
	const mfvBaseObject *hitObject = NULL;
	while (1) {
		uint32_t o = cell[2] * resolution[0] * resolution[1] + cell[1] * resolution[0] + cell[0];
		if (cells[o] != NULL) {
			cells[o]->intersect(ray, &hitObject);
			if (hitObject != NULL) { ray.color = cells[o]->color; }
		}
		uint8_t k =
			((nextCrossingT[0] < nextCrossingT[1]) << 2) +
			((nextCrossingT[0] < nextCrossingT[2]) << 1) +
			((nextCrossingT[1] < nextCrossingT[2]));
		static const uint8_t map[8] = {2, 1, 2, 1, 2, 2, 0, 0};
		uint8_t axis = map[k];

		if (ray.tmax < nextCrossingT[axis]) break;
		cell[axis] += step[axis];
		if (cell[axis] == exit[axis]) break;
		nextCrossingT[axis] += deltaT[axis];
	}

	return hitObject;
}

void Grid::render(void)
{
	//mfp_debug("Grid: render %s\n",this->name);
	// Init handle to GL object
	//if (! object) object = ViewManagerInst().GetGraphicData().GetRootObject().AddNewObject(this->name);
	// init vertex accumulator
	VtxAccumulator& vtxAccum(object.GetCurrentVtxAccumulator());
	//vtxAccum.UseVBO() = ViewManagerInst().IsVBOUsed();

	const CVec3f v(this->bbox.mn);
	float sizes[3];
	this->bbox.majorAxis(sizes);
	// Add grid nodes in ZY planes
	uint_t totalNodes(0);
	for (uint_t i = 0; i <= resolution[0]; ++i) {
		for (uint_t j = 0; j <= resolution[1]; ++j) {
			for (uint_t k = 0; k < 2; ++k) {
				CVec3f nd = v + CVec3f(cellDimension.x*i,cellDimension.y*j,sizes[2]*k);
				object.AddMeshNode(nd);
				//std::cout<< "node(" << i << ", " << j << " ," << k << ") = " << nd << std::endl;
			}
			const uint_t ed[] = {totalNodes++,totalNodes++};
			object.AddEdge(ed);
		}
	}
	// Add grid nodes in XY planes
	for (uint_t k = 0; k <= resolution[2]; ++k) {
		for (uint_t i = 0; i <= resolution[0]; ++i) {
			for (uint_t j = 0; j < 2; ++j) {
				CVec3f nd = v + CVec3f(cellDimension.x*i,sizes[1]*j,cellDimension.z*k);
				object.AddMeshNode(nd);
				//std::cout<< "node(" << i << ", " << j << " ," << k << ") = " << nd << std::endl;
			}
			const uint_t ed[] = {totalNodes++,totalNodes++};
			object.AddEdge(ed);
		}
	}
	// Add grid nodes in XZ planes
	for (uint_t j = 0; j <= resolution[1]; ++j) {
		for (uint_t k = 0; k <= resolution[2]; ++k) {
			for (uint_t i = 0; i < 2; ++i) {
				CVec3f nd = v + CVec3f(sizes[0]*i,cellDimension.y*j,cellDimension.z*k);
				object.AddMeshNode(nd);
				//std::cout<< "node(" << i << ", " << j << " ," << k << ") = " << nd << std::endl;
			}
			const uint_t ed[] = {totalNodes++,totalNodes++};
			object.AddEdge(ed);
		}
	}
	//mfp_debug("After grid render\n");
	//FV_CHECK_ERROR_GL();

}

void Grid::drawOn(bool On)
{
	//mfp_debug("DrawOn\n");
	object.GetCurrentVtxAccumulator().aRenderEdges = On;
}

void Grid::dump(void) const
{

}

/*************************************************/
int create_grid(double targetoccupancy,
		        uint32_t numelems,
		        const BBox3D& gbox,
		        const BBox3D* boxes,
		        cl_int** griddata,
		        cl_int** tridata,
		        gridinfo_t* gridinfo
		        )
{
	using namespace fvmath;
	cl_int *lgriddata, *ltridata;
	gridinfo->minDimension.x = gbox.mn.x;
	gridinfo->minDimension.y = gbox.mn.y;
	gridinfo->minDimension.z = gbox.mn.z;
	gridinfo->maxDimension.x = gbox.mx.x;
	gridinfo->maxDimension.y = gbox.mx.y;
	gridinfo->maxDimension.z = gbox.mx.z;

	float xlength = gridinfo->maxDimension.x-gridinfo->minDimension.x;
	float ylength = gridinfo->maxDimension.y-gridinfo->minDimension.y;
	float zlength = gridinfo->maxDimension.z-gridinfo->minDimension.z;
	float norm = std::cbrt((numelems/targetoccupancy)/(xlength*ylength*zlength));
	gridinfo->cellCount.x=fmax(1,std::roundf(xlength*norm));
	gridinfo->cellCount.y=fmax(1,std::roundf(ylength*norm));
	gridinfo->cellCount.z=fmax(1,std::roundf(zlength*norm));
	gridinfo->cellSize.x=xlength/gridinfo->cellCount.x;
	gridinfo->cellSize.y=ylength/gridinfo->cellCount.y;
	gridinfo->cellSize.z=zlength/gridinfo->cellCount.z;

	size_t gridsize = gridinfo->cellCount.x * gridinfo->cellCount.y * gridinfo->cellCount.z;
	//printf("cellCOunts= %d %d %d\tgridsize = %d\n",gridinfo->cellCount.x,gridinfo->cellCount.y,gridinfo->cellCount.z, gridsize);
	lgriddata = new cl_int [gridsize+1];
	memset(lgriddata,0x0,sizeof(cl_int)*(gridsize+1));

	//cl_int3 minidx;
	//cl_int3 maxidx;
	CVec3i minidx, maxidx;
	for (uint32_t i=0;i<numelems;i++)
	{
		//count into griddata
		minidx.x=fv_min(gridinfo->cellCount.x-1,(int)std::floor(((boxes[i].mn.x-gridinfo->minDimension.x)*gridinfo->cellCount.x)/xlength));
		maxidx.x=fv_min(gridinfo->cellCount.x-1,(int)std::floor(((boxes[i].mx.x-gridinfo->minDimension.x)*gridinfo->cellCount.x)/xlength));
		minidx.y=fv_min(gridinfo->cellCount.y-1,(int)std::floor(((boxes[i].mn.y-gridinfo->minDimension.y)*gridinfo->cellCount.y)/ylength));
		maxidx.y=fv_min(gridinfo->cellCount.y-1,(int)std::floor(((boxes[i].mx.y-gridinfo->minDimension.y)*gridinfo->cellCount.y)/ylength));
		minidx.z=fv_min(gridinfo->cellCount.z-1,(int)std::floor(((boxes[i].mn.z-gridinfo->minDimension.z)*gridinfo->cellCount.z)/zlength));
		maxidx.z=fv_min(gridinfo->cellCount.z-1,(int)std::floor(((boxes[i].mx.z-gridinfo->minDimension.z)*gridinfo->cellCount.z)/zlength));
		for (int x=minidx.x;x<=maxidx.x;x++)
			for (int y=minidx.y;y<=maxidx.y;y++)
				for (int z=minidx.z;z<=maxidx.z;z++)
				{
					//if (cell_triangle_intersect(triangles[i],
					//	(gridinfo->minDimension.x+(x*gridinfo->cellSize.x)),(gridinfo->minDimension.x+((x+1)*gridinfo->cellSize.x)),
					//	(gridinfo->minDimension.y+(y*gridinfo->cellSize.y)),(gridinfo->minDimension.y+((y+1)*gridinfo->cellSize.y)),
					//	(gridinfo->minDimension.z+(z*gridinfo->cellSize.z)),(gridinfo->minDimension.z+((z+1)*gridinfo->cellSize.z)) ) == 1)
						(lgriddata[x + (gridinfo->cellCount.x * y) + (gridinfo->cellCount.x * gridinfo->cellCount.y * z)])++;
				}
 	}

	for (uint32_t i=1; i<=gridsize; i++) lgriddata[i] += lgriddata[i-1];
    ltridata = new cl_int[ lgriddata[gridsize] ];

    for (uint32_t i=0;i<numelems;i++)
	{
    	//put indexes into tridata, decrease griddata
    	minidx.x=fv_min(gridinfo->cellCount.x-1,(int)std::floor(((boxes[i].mn.x-gridinfo->minDimension.x)*gridinfo->cellCount.x)/xlength));
    	maxidx.x=fv_min(gridinfo->cellCount.x-1,(int)std::floor(((boxes[i].mx.x-gridinfo->minDimension.x)*gridinfo->cellCount.x)/xlength));
    	minidx.y=fv_min(gridinfo->cellCount.y-1,(int)std::floor(((boxes[i].mn.y-gridinfo->minDimension.y)*gridinfo->cellCount.y)/ylength));
    	maxidx.y=fv_min(gridinfo->cellCount.y-1,(int)std::floor(((boxes[i].mx.y-gridinfo->minDimension.y)*gridinfo->cellCount.y)/ylength));
    	minidx.z=fv_min(gridinfo->cellCount.z-1,(int)std::floor(((boxes[i].mn.z-gridinfo->minDimension.z)*gridinfo->cellCount.z)/zlength));
    	maxidx.z=fv_min(gridinfo->cellCount.z-1,(int)std::floor(((boxes[i].mx.z-gridinfo->minDimension.z)*gridinfo->cellCount.z)/zlength));
		for (int x=minidx.x;x<=maxidx.x;x++)
			for (int y=minidx.y;y<=maxidx.y;y++)
				for (int z=minidx.z;z<=maxidx.z;z++)
				{
				//	if (cell_triangle_intersect(triangles[i],
				//			(gridinfo->minDimension.x+(x*gridinfo->cellSize.x)),(gridinfo->minDimension.x+((x+1)*gridinfo->cellSize.x)),
				//			(gridinfo->minDimension.y+(y*gridinfo->cellSize.y)),(gridinfo->minDimension.y+((y+1)*gridinfo->cellSize.y)),
				//			(gridinfo->minDimension.z+(z*gridinfo->cellSize.z)),(gridinfo->minDimension.z+((z+1)*gridinfo->cellSize.z)) ) == 1)
					lgriddata[x + (gridinfo->cellCount.x*y) + (gridinfo->cellCount.x*gridinfo->cellCount.y * z)]--;
					ltridata[lgriddata[x + (gridinfo->cellCount.x*y) + (gridinfo->cellCount.x*gridinfo->cellCount.y * z)]]=i;
				}
	}
	*griddata=lgriddata;
	*tridata=ltridata;
	return 0;
}

}
