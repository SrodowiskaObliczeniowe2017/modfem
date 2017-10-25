#define EPSILON

/* trace grid */
void trace_grid( __private ray_t local_ray, 
                 __global  hit_t * target_hit, 
		 __global const elem_info_t * triangles, 
                 __global const gridaccel_t * accel )
{
 __private gridaccel_t local_accel;
 __private gridinfo_t gridinfo;
 local_accel = *accel;
 gridinfo = *(local_accel.gridinfo);
 __private hit_t hit;
 __private float t=0,maxt=INFINITY;
 __private float3 d=fabs(gridinfo.cellSize/local_ray.direction);
 hit.elindex = -1;
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
