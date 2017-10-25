/*
 * mic_operator.cpp
 *
 *  Created on: 16 lis 2015
 *      Author: pmaciol
 */
#include "mic.h"
#include "mic_operator.h"
#include "Log.h"
#include <unistd.h>
#include <omp.h>

namespace FemViewer {
namespace MIC {

template< typename TScalar = double >
void dummy_routine(int target = mic_CPU)
{
    double *a, *b, *c;
    int i,j,k, ok;
    const int n=1000, nn=1000000;

    // allocated memory on the heap aligned to 64 byte boundary
    ok = posix_memalign((void**)&a, 64, nn*sizeof(TScalar));
    ok = posix_memalign((void**)&b, 64, nn*sizeof(TScalar));
    ok = posix_memalign((void**)&c, 64, nn*sizeof(TScalar));

    // initialize matrices
    for (i=0; i < nn; i++) a[i] = (TScalar)i;
    for (i=0; i < nn; i++) b[i] = (TScalar)i;
    for (i=0; i < nn; i++) c[i] = TScalar(0.0);

    //if (target == mic_XEONPHI) {
#ifdef XEON_PHI
          //offload code
    #pragma offload target(mi:0) in(a,b:length(nn)) inout(c:length(nn))
#endif
    	{
    	//parallelize via OpenMP on MIC
	#pragma omp parallel for
        for(i = 0; i < n; i++) {
          for(k = 0; k < n; k++) {
    #pragma vector aligned
    #pragma ivdep
        	  for( j = 0; j < n; j++ ) {
        		  //c[i][j] = c[i][j] + a[i][k]*b[k][j];
        		  c[i*n+j] = c[i*n+j] + a[i*n+k]*b[k*n+j];
        	  }
          }
        }

    	}
//    }
//    else { // code for host CPU
//
//    }


}


int MicOperator::initMIC(bool quiet)
{
	// Check the number of available cores
	long num_cores = sysconf(_SC_NPROCESSORS_ONLN);
	if (!quiet)
		mfp_log_info("Number of available cores on HOST: %d", num_cores);

	int is_xphi = 0;
#pragma offload target(mic)
	{
		long num_cores = sysconf(_SC_NPROCESSORS_ONLN);
		is_xphi++;
	}

	if (!quiet && is_xphi)
		mfp_log_info("Number of available cores on ACCELERATOR: %d", num_cores);

	// Invoke start up test procedure
	dummy_routine(mic_CPU);
	if (is_xphi) dummy_routine(mic_XEONPHI);

	return 0;
}

void MicOperator::shutdownMIC()
{
	mfp_log_debug("Shutdown procedure");
}

}
}



