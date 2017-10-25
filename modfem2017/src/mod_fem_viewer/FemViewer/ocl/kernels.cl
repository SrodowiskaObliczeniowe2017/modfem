#include "types.h"
#include "cl/grid.cl"
#include "cl/dummy.cl"

#include "types.h"
#include "cl/algs.h"

//fix a bug with float comparison
#define EPSILON 0.00001f

__kernel void finalize_grid( __global struct accel_s * accelerator, __global const struct gridinfo_s * gridinfo,
			__global const int * grid, __global const int * cells)
{
 accelerator->type=ACCEL_GRID;
 accelerator->grid.gridinfo = gridinfo;
 accelerator->grid.grid = grid;
 accelerator->grid.cells = cells;
}

void trace_grid( __private struct ray_s local_ray, __global struct hit_s * target_hit,
			__global const struct triangle_s * triangles, __global const struct gridaccel_s * accel )
{
 __private struct gridaccel_s local_accel;
 __private struct gridinfo_s gridinfo;
 local_accel = *accel;
 gridinfo = *(local_accel.gridinfo);
 __private struct hit_s hit;
 __private float t=0,maxt=INFINITY;
 __private float3 d=fabs(gridinfo.cellSize/local_ray.direction);
 hit.triIndex= -1;
 hit.rayDirection=normalize(local_ray.direction);
 hit.reconstructInfo=local_ray.reconstructInfo;
 target_hit->reconstructInfo=local_ray.reconstructInfo;



 //put ray origin into the grid
 __private float mintx=(local_ray.direction.x>0.0f)?
	(gridinfo.minDimension.x-local_ray.origin.x)/local_ray.direction.x:
	(gridinfo.maxDimension.x-local_ray.origin.x)/local_ray.direction.x;
 __private float minty=(local_ray.direction.y>0.0f)?
	(gridinfo.minDimension.y-local_ray.origin.y)/local_ray.direction.y:
	(gridinfo.maxDimension.y-local_ray.origin.y)/local_ray.direction.y;
 __private float mintz=(local_ray.direction.z>0.0f)?
	(gridinfo.minDimension.z-local_ray.origin.z)/local_ray.direction.z:
	(gridinfo.maxDimension.z-local_ray.origin.z)/local_ray.direction.z;
 t=fmax(0.0f,fmax(mintx,fmax(minty,mintz)));
 local_ray.origin+=(t+EPSILON)*local_ray.direction;
// !!!BIG FAT FLOATING POINT ISSUE - origin not in grid after correction!!!
 t=EPSILON;

 //set maxt to grid exit point
 __private float maxtx=(local_ray.direction.x>0)?
	(gridinfo.maxDimension.x-local_ray.origin.x)/local_ray.direction.x:
	(gridinfo.minDimension.x-local_ray.origin.x)/local_ray.direction.x;
 __private float maxty=(local_ray.direction.y>0)?
	(gridinfo.maxDimension.y-local_ray.origin.y)/local_ray.direction.y:
	(gridinfo.minDimension.y-local_ray.origin.y)/local_ray.direction.y;
 __private float maxtz=(local_ray.direction.z>0)?
	(gridinfo.maxDimension.z-local_ray.origin.z)/local_ray.direction.z:
	(gridinfo.minDimension.z-local_ray.origin.z)/local_ray.direction.z;
 maxt=fmin(maxtx,fmin(maxty,maxtz));

 if (maxt<t)
 {
  //misses the whole grid
  *target_hit=hit;
  return;
 }


 //find the cell containing ray origin
 __private int3 current=min(gridinfo.cellCount-1,convert_int3(floor(
	((local_ray.origin-gridinfo.minDimension)/(gridinfo.maxDimension-gridinfo.minDimension))
	*convert_float3(gridinfo.cellCount))));

 __private int3 pitch;

 pitch.x = 1;
 pitch.y = gridinfo.cellCount.x;
 pitch.z = gridinfo.cellCount.x*gridinfo.cellCount.y;



// __private float3 dt=remainder((local_ray.origin-gridinfo.minDimension),gridinfo.cellSize)/(-local_ray.direction);
 __private float3 dt=fmod((local_ray.origin-gridinfo.minDimension),gridinfo.cellSize)/(-local_ray.direction);
 if (dt.x<0.0f)
  dt.x+=d.x;
 if (dt.y<0.0f)
  dt.y+=d.y;
 if (dt.z<0.0f)
  dt.z+=d.z;
 __private float tt;
 __private float3 tbary;
 __private int cellindex;
 __private int currenttri;
 __private int nextcellindex;
 __private struct triangle_s ttri;

 //traverse the grid until we have a hit
 while (maxt-t>EPSILON)
 {
//  cellindex = dot(current,pitch);
  cellindex = (current.x*pitch.x + current.y*pitch.y + current.z*pitch.z);
if( (current.x<0)||(current.y<0)||(current.z<0)||(current.x>=gridinfo.cellCount.x)||(current.y>=gridinfo.cellCount.y)||(current.z>=gridinfo.cellCount.z))
{
hit.triIndex=-2;
hit.rayDirection=dt;
hit.baryCoords=local_ray.origin;
break;
}
  //manual cacheing in private memory
  nextcellindex = local_accel.grid[cellindex+1];
  for(__private int i = local_accel.grid[cellindex];i<nextcellindex;i++)
  {
   ttri=triangles[currenttri=local_accel.cells[i]];
   tt=intersectRayTri(&local_ray,&ttri,&tbary);
   if ((tt<maxt)&&(tt>t))
   {
    maxt=tt;
    hit.triIndex=currenttri;
    hit.baryCoords=tbary;
   }
  }
  //step into next cell
  __private float dmin=fmin(dt.x,fmin(dt.y,dt.z));
  t+=dmin;
  dt-=dmin;
  if(dt.x<=0.0f)
  {
   dt.x+=d.x;
   if (local_ray.direction.x>0.0f)
    current.x++;
   else
    current.x--;
  }
  if(dt.y<=0.0f)
  {
   dt.y+=d.y;
   if (local_ray.direction.y>0.0f)
    current.y++;
   else
    current.y--;
  }
  if(dt.z<=0.0f)
  {
   dt.z+=d.z;
   if (local_ray.direction.z>0.0f)
    current.z++;
   else
    current.z--;
  }
 }

 *target_hit=hit;
 return;
}



bool trace_shadow_grid(__private float3 origin, __private float3 direction, __private float maxt,
        __global const struct triangle_s * triangles, __global const struct gridaccel_s * accel )
{
 __private struct gridaccel_s local_accel;
 __private struct gridinfo_s gridinfo;
 __private struct ray_s local_ray;
 local_ray.origin = origin;
 local_ray.direction = direction;
 local_accel = *accel;
 gridinfo = *(local_accel.gridinfo);
 __private float t=0;
 __private float3 d=fabs(gridinfo.cellSize/local_ray.direction);

 //put ray origin into the grid
 __private float mintx=(local_ray.direction.x>0.0f)?
        (gridinfo.minDimension.x-local_ray.origin.x)/local_ray.direction.x:
        (gridinfo.maxDimension.x-local_ray.origin.x)/local_ray.direction.x;
 __private float minty=(local_ray.direction.y>0.0f)?
        (gridinfo.minDimension.y-local_ray.origin.y)/local_ray.direction.y:
        (gridinfo.maxDimension.y-local_ray.origin.y)/local_ray.direction.y;
 __private float mintz=(local_ray.direction.z>0.0f)?
        (gridinfo.minDimension.z-local_ray.origin.z)/local_ray.direction.z:
        (gridinfo.maxDimension.z-local_ray.origin.z)/local_ray.direction.z;
 t=fmax(0.0f,fmax(mintx,fmax(minty,mintz)));
 local_ray.origin+=(t+EPSILON)*local_ray.direction;
// !!!BIG FAT FLOATING POINT ISSUE - origin not in grid after correction!!!
 t=EPSILON;

 //set maxt to grid exit point
 __private float maxtx=(local_ray.direction.x>0)?
        (gridinfo.maxDimension.x-local_ray.origin.x)/local_ray.direction.x:
        (gridinfo.minDimension.x-local_ray.origin.x)/local_ray.direction.x;
 __private float maxty=(local_ray.direction.y>0)?
        (gridinfo.maxDimension.y-local_ray.origin.y)/local_ray.direction.y:
        (gridinfo.minDimension.y-local_ray.origin.y)/local_ray.direction.y;
 __private float maxtz=(local_ray.direction.z>0)?
        (gridinfo.maxDimension.z-local_ray.origin.z)/local_ray.direction.z:
        (gridinfo.minDimension.z-local_ray.origin.z)/local_ray.direction.z;
 maxt=fmin(maxt,fmin(maxtx,fmin(maxty,maxtz)));

 if (maxt<t)
 {
  return true;
 }

 //find the cell containing ray origin
 __private int3 current=min(gridinfo.cellCount-1,convert_int3(floor(
        ((local_ray.origin-gridinfo.minDimension)/(gridinfo.maxDimension-gridinfo.minDimension))
        *convert_float3(gridinfo.cellCount))));

 __private int3 pitch;
 pitch.x = 1;
 pitch.y = gridinfo.cellCount.x;
 pitch.z = gridinfo.cellCount.x*gridinfo.cellCount.y;

 __private float3 dt=fmod((local_ray.origin-gridinfo.minDimension),gridinfo.cellSize)/(-local_ray.direction);
 if (dt.x<0.0f)
  dt.x+=d.x;
 if (dt.y<0.0f)
  dt.y+=d.y;
 if (dt.z<0.0f)
  dt.z+=d.z;
 __private float tt;
 __private float3 tbary;
 __private int cellindex;
 __private int nextcellindex;
 __private struct triangle_s ttri;

 while (maxt-t>EPSILON)
 {
//  cellindex = dot(current,pitch);
  cellindex = (current.x*pitch.x + current.y*pitch.y + current.z*pitch.z);
if( (current.x<0)||(current.y<0)||(current.z<0)||(current.x>=gridinfo.cellCount.x)||(current.y>=gridinfo.cellCount.y)||(current.z>=gridinfo.cellCount.z))
break;
  nextcellindex = local_accel.grid[cellindex+1];
  for(__private int i = local_accel.grid[cellindex];i<nextcellindex;i++)
  {
   ttri=triangles[local_accel.cells[i]];
   tt=intersectRayTri(&local_ray,&ttri,&tbary);
   if ((tt<maxt)&&(tt>t))
   {
    return false;
   }
  }
  //step into next cell
  __private float dmin=fmin(dt.x,fmin(dt.y,dt.z));
  t+=dmin;
  dt-=dmin;
  if(dt.x<=0.0f)
  {
   dt.x+=d.x;
   if (local_ray.direction.x>0.0f)
    current.x++;
   else
    current.x--;
  }
  if(dt.y<=0.0f)
  {
   dt.y+=d.y;
   if (local_ray.direction.y>0.0f)
    current.y++;
   else
    current.y--;
  }
  if(dt.z<=0.0f)
  {
   dt.z+=d.z;
   if (local_ray.direction.z>0.0f)
    current.z++;
   else
    current.z--;
  }
 }

 return true;
}


void trace(__private struct ray_s ray, __global struct hit_s * target_hit,
	__global const struct triangle_s * triangles, __global const struct accel_s * accel)
{

 switch (accel->type)
 {
  case ACCEL_DUMMY:
   trace_dummy(ray, target_hit, triangles, (__global struct dummyaccel_s *) &(accel->dummy));
   return;
  case ACCEL_GRID:
   trace_grid(ray, target_hit, triangles, (__global struct gridaccel_s *) &(accel->grid));
   return;
  default:
   return;
 }

}

bool trace_shadow(__private float3 origin, __private float3 direction, __private float maxt,
	__global const struct triangle_s * triangles, __global const struct accel_s * accel)
{
 switch (accel->type)
 {
  case ACCEL_DUMMY:
   return trace_shadow_dummy(origin, direction, maxt, triangles, (__global struct dummyaccel_s *) &(accel->dummy));
  case ACCEL_GRID:
   return trace_shadow_grid(origin, direction, maxt, triangles, (__global struct gridaccel_s *) &(accel->grid));
 }
 return false;
}

/* kernel for tracing rays*/
__kernel void kernel_cast_ray(
		__global const struct ray_s * rays,
		int raycount, __global struct hit_s * hits,
		__global const struct triangle_s * triangles,
		__global const struct accel_s * accel)
{
	__private int raynum = get_global_id(0)+get_global_id(1)*get_global_size(0);
	if (raynum >= raycount)
		return;
	trace(rays[raynum], &(hits[raynum]),triangles, accel);

	return;
}
