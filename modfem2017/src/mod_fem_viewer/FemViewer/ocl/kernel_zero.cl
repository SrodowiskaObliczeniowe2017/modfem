__kernel void zero( __global float *vh)
{
  const float4 id = (float4)(get_global_id(0), get_global_id(1), get_global_id(2), 0);
  const float4 sz = (float4)(1, get_global_size(0), get_gloabal_size(0)*get_gloabal_size(1), 0);

  vh[(int)dot(id,sz)] = 0;
}

// Matrix multiplication kernel called by MatrixMul()
__kernel void MatVecMul( const __global float* M,
                         const __global float* V,
                         uint width, uint height,
                         __global float* W)
{
    // Row index
    uint y = get_global_id(0);
    if (y < height) {
    
        // Row pointer
        const __global float* row = M + y * width;

        // Compute dot product  
        float dotProduct = 0;
        for (uint x = 0; x < width; ++x)
            dotProduct += row[x] * V[x];

        // Write result to global memory
        W[y] = dotProduct;
    }
}
