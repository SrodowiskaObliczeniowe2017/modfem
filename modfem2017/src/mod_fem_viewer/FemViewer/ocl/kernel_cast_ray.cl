
struct julia_context {
  float4 dir_top_start;
  float4 dir_bottom_start;
  float4 dir_bottom_stop;
  float4 eyeP;
  float4 lightP;
  int2 window_size;
  float epsilon;
  int maxIterations;
  int stride;
  int pad[3];
};

/* kernel for casting rays*/
__kernel void kernel_cast_ray(
	__global const ray_t * rays,
        int raycount, 
        __global hit_t * hits,
	__global const elem_info * elems,
	__global const coeffs_info * accel)
{
	__private int raynum = get_global_id(0)+get_global_id(1)*get_global_size(0);

	if (raynum >= raycount)
		return;
	
	trace(rays[raynum], &(hits[raynum]),triangles, accel);
	trace_grid(rays[raynum], &(hits[raynum]), elems, (__global gridaccel_s *) &(accel->grid));

	return;
}


__kernel void gen_primary_rays (__global struct ray_s * rays, float3 eyepoint, float3 lookat, float3 upvector, float2 angle, int2 curtile, int2 maxtile, __global unsigned int * seeds)
{
 __private int raynum = get_global_id(0) + get_global_id(1)*get_global_size(0);
 __local float3 imageplanenormal;
 imageplanenormal = eyepoint-lookat;
 __local float3 xdest;
 __local float3 ydest;
 __local float2 targetlength;
 xdest = cross(upvector,imageplanenormal);
 ydest = cross(imageplanenormal,xdest);
 targetlength.x = tan(deg2rad(angle.x)/2)*length(imageplanenormal);
 targetlength.y = tan(deg2rad(angle.y)/2)*length(imageplanenormal);
 ydest=ydest*(targetlength.y/length(ydest));
 xdest=xdest*(targetlength.x/length(xdest));
 __private struct ray_s local_ray;
 local_ray.origin=eyepoint;
 local_ray.direction= -imageplanenormal;
 local_ray.direction+=(((( (get_global_id(0)+curtile.x*get_global_size(0)) +0.5f)/(get_global_size(0)*maxtile.x))*2)-1.0f)*xdest;
 local_ray.direction+=(((( (get_global_id(1)+curtile.y*get_global_size(1)) +0.5f)/(get_global_size(1)*maxtile.y))*2)-1.0f)*ydest;
 local_ray.direction=normalize(local_ray.direction);
 local_ray.reconstructInfo.randomSeed=seeds[raynum];
 local_ray.reconstructInfo.pixelNumber=raynum;
 local_ray.reconstructInfo.importance=1.0f;
 rays[raynum]=local_ray;
}

__kernel void copy_image_to_pbo(__global const struct pixel_s * image, int2 curtile, int rowsize, __global uchar4 * pbo)
{
 __private int incoord = get_global_id(1)*get_global_size(0) + get_global_id(0);
}
