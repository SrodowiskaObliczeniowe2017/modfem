#include"ExponentEngine.hpp"

//static ExponentEngine<int,double> engine;

void multiply_1d(int pdeg,double pt,double* outPhi,double* outdPhi);

int shape_functions_3D(int pdeg,int base,double* pt,
		double* Phi,double* dPhix,double* dPhiy,double* dPhiz)
{
	//static ExponentEngine<int,double> eg;
	double phi[3][10], dphi[3][10];
	int pdeg0,pdeg1,pdeg2;
	pdeg0 = pdeg / 100; pdeg1 = pdeg2 = pdeg0;
	//int pdeg0 = degree_3d(pdeg,base,0);
	multiply_1d(pdeg0,pt[0],phi[0],dphi[0]);
	pdeg0++;
	//int pdeg1 = degree_3d(pdeg,base,1);
	multiply_1d(pdeg1,pt[1],phi[1],dphi[1]);
	pdeg1++;
	//int pdeg2 = degree_3d(pdeg,base,2);
	multiply_1d(pdeg2,pt[2],phi[2],dphi[2]);
	pdeg2++;
	//printf("pdeg: %d %d %d\n",pdeg0,pdeg1,pdeg2);
	int kk = 0, i,j,k;
	if (base == eTensor)
	{

		for( i = 0; i < pdeg0; ++i) //phiz
		{
			double z = phi[2][i];
			for ( j = 0; j < pdeg1-i; ++j) //phiy
			{
				double y = z*phi[1][j];
				for( k = 0; k < pdeg2-j-i; ++k) //phix
	            {
					Phi[kk++] = y*phi[0][k];
	            }
	        }
		}
		if (dPhix != NULL) {
			kk = 0;
			for( i = 0; i < pdeg0; ++i) //phiz
			{
				double z = phi[2][i];
				double dz = dphi[2][i];
				for ( j = 0; j < pdeg1-i; ++j) //phiy
				{
					double y = phi[1][j];
					double yz = y * z;
					double dy = dphi[1][j];
					for( k = 0; k < pdeg2-j-i; ++k) //phix
					{
						double x = phi[0][k];
						double dx = dphi[0][k];
						dPhix[kk] = dx*yz;
						dPhiy[kk] = dy*x*z;
						dPhiz[kk++] = dz*x*y;
					}
				}
			}

		}
	}
	else { // eTensor
		//printf("dla tensor\n");
		for ( i = 0; i < pdeg0; ++i) //x
	    {
			double x = phi[0][i];
	        for ( j = 0; j < pdeg1; ++j) //z
	        {
	        	double z = x*phi[2][j];
	           	for( k = 0; k < pdeg2-i; ++k) //y
	           	{
	           		Phi[kk++] = z*phi[1][k];
	            }
	        }
		}
		if (dPhix != NULL) {
			kk = 0;
			for( i = 0; i < pdeg0; ++i) //phix
			{
				double x = phi[0][i];
				double dx = dphi[0][i];
				for ( j = 0; j < pdeg1; ++j) //phiz
				{
					double z = phi[2][j];
					double zx = z*x;
					double dzx = x*dphi[2][j];
					double dxz = dx * z;
					for( k = 0; k < pdeg2-i; ++k) //phiy
					{
						double y = phi[1][k];
						double dy = dphi[1][k];
						dPhix[kk] = dxz*y;
						dPhiy[kk] = dy*zx;
						dPhiz[kk++] = dzx*y;
					}
				}
			}
		}
	}// else eTensor
	//printf("num_shap_kk = %d\n",kk);
	return kk;
}


int degree_3d(int pdeg,int base_type,int direction)
{
	int pdegx,pdegz;
	if (base_type == eTensor) {
       		pdegx = pdeg % 100;
       		pdegz = pdeg / 100;
		if (direction == 0) return (pdegx % 10);
		if (direction == 1) return (pdegx % 10);
		return pdegz;
     	}
     	else if (base_type == eComplete) {
       		if (pdeg < 100) {
	     		pdegx = pdeg;
 	     		pdegz = pdeg;
			if (direction == 2) return pdegz;
			return pdegx;
       		}
       		else {
	     		pdegx = pdeg % 100;
	     		pdegz = pdeg / 100;
			if (direction == 0) return (pdegx % 10);
			if (direction == 1) return (pdegx % 10);
			return pdegz;
       		}
     	}
     	else {
       		printf("Type of base not valid for prisms!\n");
     	}
	return 0;
}

void multiply_1d(int pdeg,
		double pt,double* outPhi,double* outdPhi)
{
	outPhi[0] = 1.0;
	if (pdeg > 0)  outPhi[1] = pt;
	int i;
	for ( i=2;i<=pdeg;++i) {
		outPhi[i] = pt*outPhi[i-1];
	}

	if (outdPhi != NULL) {
		outdPhi[0] = 0.0;
		if (pdeg > 0) outdPhi[1] = 1.0;
		for ( i=2;i<=pdeg;++i) {
			outdPhi[i] = pt * outdPhi[i-1];
		}
		for ( i=2;i<=pdeg;++i) {
			outdPhi[i] *= i;
		}
	}
}

