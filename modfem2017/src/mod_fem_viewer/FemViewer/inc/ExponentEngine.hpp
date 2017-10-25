#ifndef __Exponent_Engine_hpp__
#define __Exponent_Engine_hpp__

#include<stdio.h>
#include<stdlib.h>
#include<assert.h>
#include<algorithm>
#include<iostream>
#include<stdexcept>
#include<vector>
//#include"aph_intf.h"

#define fv_max(a,b) (a) > (b) ? (a) : (b)

enum eBaseType {
	eTensor = 0x0,//APC_BASE_TENSOR_DG,
	eComplete = 0x1,//APC_BASE_COMPLETE_DG,
};


extern int shape_functions_3D(int pdeg,int base,double*,double*,double*,double*,double*);
int degree_3d(int pdeg,int base_type,int direction);
template<typename TExp,typename TScalar> class ExponentEngine;

template<typename TExp,typename TScalar>
std::ostream& operator << (std::ostream& os, const ExponentEngine<TExp,TScalar>& rhs);

template<typename TExp,typename TScalar>
class ExponentEngine 
{
  public:
	struct fvType2 {
		union {
			TExp u[2];
			struct { TExp x,y; };
		};
	};
	struct fvType3 {
		union {
			TExp u[3];
			struct { TExp x,y,z; };
		};
		bool operator==(const fvType3& rhs) {
			return (x == rhs.x && y == rhs.y && z == rhs.z);
		}
	};

	typedef std::vector<fvType3> vExpm;

  public:
	ExponentEngine(int pdeg = 101,int base = 0);
	~ExponentEngine();
	
	void Set(int pdeg,int base);
	const fvType3* GetExpMap() const { return _expms.data(); }
	int GetPdeg() const { return _pdeg; }
	int GetBaseType() const { return _base; }
	int GetNumShapeFunc() const { return _num_shap; }

  protected:
	vExpm     _expms;
	int 	  _pdeg;
	eBaseType _base;
	int       _num_shap;

  private:
	void init();
        bool test();
	int  convert(vExpm& expms);
	static int degree_3d(int pdeg,int base_type, int direction);
	static int num_shape_func(int pdeg,int base);
	static int order_3d(int pdeg,int base_type, fvType3 expm[]);
	static int order_2d(int pdeg,int base_type, fvType2 expm[]);
    static int order_1d(int pdeg,TExp expm[]);
    static void multiply_1d(int pdeg,TScalar pt,TScalar* outPhi,TScalar* outdPhi);

    ExponentEngine(const ExponentEngine&);
    ExponentEngine& operator=(const ExponentEngine&);

    //friend std::ostream& operator<< <>(std::ostream& os,const ExponentEngine& e);
    friend int shape_functions_3D(int pdeg,int base,double*,double*,double*,double*,double*);
};

template<typename TExp,typename TScalar>
ExponentEngine<TExp,TScalar>::ExponentEngine(int pdeg,int base) {
	Set(pdeg,base);
}

template<typename TExp,typename TScalar>
ExponentEngine<TExp,TScalar>::~ExponentEngine() {
	_expms.clear();
}


template<typename TExp,typename TScalar>
void ExponentEngine<TExp,TScalar>::Set(int pdeg,int base) {
	if (pdeg > 10 && base == eComplete) {
		char msg[64];
		sprintf(msg,"Warrning: degre value %d is out of range (<10) for COMPLETE type!\n",pdeg);
		throw std::runtime_error(msg);
	}
	if (pdeg < 101 && base == eTensor) {
		char msg[64];
		sprintf(msg,"Warrning: degre value %d is invalid (>=101) for TENSOR type!\n",pdeg);
		throw std::runtime_error(msg);
	}
	if (base != eTensor && base != eComplete) {
		char msg[64];
		sprintf(msg,"Warrning: unknown base type: %d allowed %d - TENSOR and %d - COMPLETE!\n",base,eTensor,eComplete);
		throw std::runtime_error(msg);
	}
	_pdeg = pdeg;
	_base = static_cast<eBaseType>(base);
	init();
	if (!test()) throw std::runtime_error("test failed!");
	vExpm expms;
	convert(expms);
	expms.clear();
}

template<typename TExp,typename TScalar>
void ExponentEngine<TExp,TScalar>::init() {
	_num_shap = num_shape_func(_pdeg,_base);
	if (_num_shap) {
		if (!_expms.empty())_expms.clear();
		_expms.resize(_num_shap);
		int num = order_3d(_pdeg,_base,_expms.data());
		printf("num = %d shape = %d\n",num,_num_shap);
		assert(num == _num_shap);
	}
}

template<typename TExp,typename TScalar>
bool ExponentEngine<TExp,TScalar>::test() {
	if (_base == eTensor) return true;
	bool passed = true;
	for (int i=0; i < _num_shap; ++i) {
		const fvType3* p = &_expms[i];
		TExp sum = p->x + p->y + p->z;
		if (sum > _pdeg) {
			passed = false;
			std::cout << "Test fsild for phi[" << i << "] = { "
			<< p->x << ", " << p->y << ", " << p->z << "}\n";
		}
	}
	return passed;
}


template<typename TExp,typename TScalar>
int ExponentEngine<TExp,TScalar>::convert(vExpm& expms)
{
	if (!expms.empty()) expms.clear();
	int pdeg0 = degree_3d(_pdeg,_base,0); pdeg0++;
	int pdeg1 = degree_3d(_pdeg,_base,1); pdeg1++;
	int pdeg2 = degree_3d(_pdeg,_base,2); pdeg2++;
	printf("pdeg: %d %d %d\n",pdeg0,pdeg1,pdeg2);
	int kk = 0;
	if (_base == eComplete) {
	for(int i = 0; i < pdeg0; ++i) //phiz
        {
                for(int j = 0; j < pdeg1-i; ++j) //phiy
                    for(int k = 0; k < pdeg2-j-i; ++k) //phix
                    {
			fvType3 el = {k,j,i};
			expms.push_back(el);
			++kk;
                     	//expm[kk].x = i;expm[kk].y = k;expm[kk++].z = j;
                    }
        }
	}
	else { // eTensor
	for(int i = 0; i < pdeg0; ++i) //phix
        {
                for(int j = 0; j < pdeg1; ++j) //phiz
                    for(int k = 0; k < pdeg2-i; ++k) //phiy
                    {
			fvType3 el = {i,k,j};
			expms.push_back(el);
			++kk;
                     	//expm[kk].x = i;expm[kk].y = k;expm[kk++].z = j;
                    }
        }
	}

	std::vector<int> pos;
	int i = 0;

		for (auto el : expms) {
			printf("el: {%d, %d, %d}\n",el.x,el.y,el.z);
			auto it = std::find(_expms.begin(),_expms.end(),el);
			int p = std::distance(_expms.begin(),it);
			//printf("p= %d\tcurr = {%d, %d, %d}\n",p,_expms[p].x,_expms[p].y,_expms[p].z);
			assert(p < _num_shap);
			pos.push_back(p);
		}



	printf("kk = %d %d\n",kk,pos.size());
assert(pos.size() == expms.size());
	for (int i =0; i < _num_shap; ++i)
		printf("%d ref = {%d, %d %d}\t%d\tcurr = {%d, %d, %d}\n",i,expms[i].x,expms[i].y,expms[i].z,pos[i],
		_expms[i].x,_expms[i].y,_expms[i].z);

	return static_cast<int>(pos.size());
}




template<typename TExp,typename TScalar>
int ExponentEngine<TExp,TScalar>::num_shape_func(int pdeg,int base) {
	int num_shap = 0;
printf("in num_shape_func pdeg = %d\n",pdeg);
    switch (base) {
    	case eTensor: {
    		int porderz,porder;
    		porderz = pdeg / 100;
    		porder = pdeg % 100;
    		num_shap = (porderz+1)*(porder+1)*(porder+2)/2;
printf("in tensor num_shape_func = %d\n",num_shap);
    	} break;
    	case eComplete: {
    		int porder=pdeg;
    		num_shap = (porder+1)*(porder+2)*(porder+3)/6;
printf("in complate num_shape_func = %d\n",num_shap);
    	} break;
    	default:
    		std::cerr << "Warrning: invalid base type for number of shape function evaluation: " << base << std::endl;
    		break;
   }
   return (num_shap);
}




template<typename TExp,typename TScalar>
int ExponentEngine<TExp,TScalar>::order_3d(int pdeg,int base_type, fvType3 expm[])
{
printf("in order_3d\n");
	if (pdeg == 0) {
	  expm[0].x = expm[1].y = expm[2].z = 0;
      return 1;
	}
	else if (pdeg == 100) {
      expm[0].x = 0; expm[0].y = 0; expm[0].z = 0;
      expm[1].x = 0; expm[1].y = 0; expm[1].z = 1;
      return 2;
    }
    else if (pdeg == 1 && base_type == eTensor) {
     expm[0].x = 0; expm[0].y = 0; expm[0].z = 0;
     expm[1].x = 1; expm[1].y = 0; expm[1].z = 0;
     expm[2].x = 0; expm[2].y = 1; expm[2].z = 0;
     return 3;   
   }
   else if (pdeg == 2 && base_type == eTensor) {
     expm[0].x = 0; expm[0].y = 0; expm[0].z = 0;
     expm[1].x = 1; expm[1].y = 0; expm[1].z = 0;
     expm[2].x = 0; expm[2].y = 1; expm[2].z = 0;
     expm[3].x = 2; expm[3].y = 0; expm[3].z = 0;
     expm[4].x = 1; expm[4].y = 1; expm[4].z = 0;
     expm[5].x = 0; expm[5].y = 2; expm[5].z = 0;
     return 6;
   }
   else {
printf("in linear and highorder_3d\n");
     expm[0].x = 0; expm[0].y = 0; expm[0].z = 0;
     expm[1].x = 1; expm[1].y = 0; expm[1].z = 0;
     expm[2].x = 0; expm[2].y = 1; expm[2].z = 0;
     expm[3].x = 0; expm[3].y = 0; expm[3].z = 1;
     if (pdeg == 1) return 4;
  
     expm[4].x = 1; expm[4].y = 0; expm[4].z = 1;
     expm[5].x = 0; expm[5].y = 1; expm[5].z = 1;
     if (pdeg == 101) return 6;
     
     expm[6].x = 1; expm[6].y = 1; expm[6].z = 0;
     if (pdeg == 2 && base_type == eComplete) {
       expm[7].x = 2; expm[7].y = 0; expm[7].z = 0;
       expm[8].x = 0; expm[8].y = 2; expm[8].z = 0;
       expm[9].x = 0; expm[9].y = 0; expm[9].z = 2;
       return 10; 
     }
     
     expm[7].x = 1; expm[7].y = 1; expm[7].z = 1;
     int base_type_2D = eComplete;
     int pdegx,pdegz;
     if (base_type == eTensor) {
       pdegx = pdeg % 100;
       pdegz = pdeg / 100;
     }
     else if (base_type == eComplete) {
       if (pdeg < 100) {
	     pdegx = pdeg;
 	     pdegz = pdeg;
       }
       else {
	     pdegx = pdeg % 100;
	     pdegz = pdeg / 100;
	 std::cout << "pdegs: " << pdegx << " " << pdegz << std::endl;
       }
     }
     else {
       printf("Type of base not valid for prisms!\n");
       exit(-1);
     }

     fvType2 ssn[60];
     TExp ttn[10];
     int num_shap_2D = order_2d(pdegx,base_type_2D,ssn);
     assert(num_shap_2D < 60);
     order_1d(pdegz,ttn);

     int num_shap;
     if (base_type == eTensor) {
       int kk = 6;
printf("in etensor\n");
       for(int j=0;j<2;j++) {
    	 for(int i=3;i<num_shap_2D;i++){
    	   expm[kk].x = ssn[i].x; expm[kk].y = ssn[i].y; expm[kk++].z = ttn[j];
         }
       }
       for(int j=2;j<pdegz+1;j++) {
         for(int i=0;i<num_shap_2D;i++) {
    	   expm[kk].x = ssn[i].x; expm[kk].y = ssn[i].y; expm[kk++].z = ttn[j];
    	 }
       }
       num_shap = kk;
     }
     else if (base_type == eComplete) {
	printf("in complate\n");
       int kk = 8, ii= 0;
       for (int j=0;j<2;j++) {
    	 for (int i=0;i<2;i++) {
    	   for (int k=2;k<=pdeg-i-j;k++){
    	      expm[kk].x = ssn[ii].x; expm[kk].y = ssn[ii].y; expm[kk++].z = ttn[k];
    	   }
    	   ii++;
    	 }
       }
       for (int j=0;j<2;j++) {
         for (int i=2;i<=pdeg-j;i++) {
    	   for(int k=0;k<=pdeg-i-j;k++){
    	     expm[kk].x = ssn[ii].x; expm[kk].y = ssn[ii].y; expm[kk++].z = ttn[k];
    	   }
    	   ii++;
    	 }
       }
       for (int j=2;j<=pdeg;j++) {
    	 for (int i=0;i<=pdeg-j;i++) {
    	   for(int k=0;k<=pdeg-i-j;k++){
    	     expm[kk].x = ssn[ii].x; expm[kk].y = ssn[ii].y; expm[kk++].z = ttn[k];
    	   }
    	   ii++;
    	 }
       }
       num_shap=kk;
     }
     else {
       printf("Unsupported base type: %d\n",base_type);
     }
     return num_shap;
   } /* end for linear and high order aproksimation */
printf("This is not correct\n");
   return(0);
}

template<typename TExp,typename TScalar>
int ExponentEngine<TExp,TScalar>::order_2d(int pdeg,int base_type,fvType2 expm[])
{
	if (pdeg == 0) {
		expm[0].x = 0; expm[0].y = 0;
		return(1);
	}
	else if (pdeg == 10) {
		expm[0].x = 0; expm[0].y = 0;
		expm[1].x = 0; expm[1].y = 1;
		return 2;
	}
	else if (pdeg == 1 && base_type == eTensor) {
		/* for linear in x and piecewise constant in y approximation */
		expm[0].x = 0; expm[0].y = 0;
		expm[1].x = 1; expm[1].y = 0;
		return 2;
	}
	else {
		/* for linear and higher order approximation */
		expm[0].x = 0; expm[0].y = 0;
		expm[1].x = 1; expm[1].y = 0;
		expm[2].x = 0; expm[2].y = 1;
		if (pdeg == 1){
			/* for linear approximation */
			return 3;
		}
		expm[3].x = 1; expm[3].y = 1;
		if (pdeg == 11) {
			/* for bi-linear approximation */
			return 4;
		}

		/* take into account different pdeg for anisotropic quads */
		int pdego = pdeg;
		int pdegx,pdegy;
		if (pdeg > 10){
			pdegx = pdeg % 10;
			pdegy = pdeg / 10;
			pdeg = fv_max(pdegx,pdegy);
		}
		else {
			pdegx = pdeg;
			pdegy = pdeg;
		}

		/* get 1D shape functions in directions x and y */
		TExp ttn[10], ssn[10];
		order_1d(pdegx,ssn);
		order_1d(pdegy,ttn);
		int ii = 4;
		if (pdego < 10) {
			if (base_type == eComplete) {
				for (int j = 0; j < 2; j++) {
					for (int i = 2; i <= pdeg-j; i++) {
						expm[ii].x = ssn[i]; expm[ii++].y = ttn[j];
					}
				}
				for (int j=2; j<= pdeg; j++) {
					for (int i=0; i<=pdeg-j; i++) {
						expm[ii].x = ssn[i]; expm[ii++].y = ttn[j];
					}
				}
			}
			else if (base_type == eTensor) {
				ii=4;
				for (int j=0;j<2;j++) {
					for (int i=2;i<=pdeg;i++) {
						expm[ii].x = ssn[i]; expm[ii++].y = ttn[j];
					}
				}
				for (int j=2;j<=pdeg;j++) {
					for (int i=0;i<=pdeg;i++) {
						expm[ii].x = ssn[i]; expm[ii++].y = ttn[j];
					}
				}
			}
		} else {
			for (int j=0;j<2;j++) {
				for (int i=2;i<=pdegx;i++) {
					expm[ii].x = ssn[i]; expm[ii++].y = ttn[j];
				}
			}
			for (int j=2;j<=pdegy;j++) {
				for (int i=0;i<=pdegx;i++) {
					expm[ii].x = ssn[i]; expm[ii++].y = ttn[j];
				}
			}
		}
		return(ii);
	} /* end if linear and higher order approximation */

	printf("Wrong parameters for order_2d!\n");
    return 0;
}

template<typename TExp,typename TScalar>
inline int ExponentEngine<TExp,TScalar>::order_1d(int pdeg,TExp expm[]) {
	/* monomial basis */
	for (int i = 0; i <= pdeg; ++i) expm[i] = i;
	return(pdeg+1);
}

template<typename TExp,typename TScalar>
inline void ExponentEngine<TExp,TScalar>::multiply_1d(int pdeg,
		TScalar pt,TScalar* outPhi,TScalar* outdPhi)
{
	outPhi[0] = 1.0;
	if (pdeg > 0)  outPhi[1] = pt;
	for (int i=2;i<=pdeg;++i) {
		outPhi[i] = pt*outPhi[i-1];
	}

	if (outdPhi != nullptr) {
		outdPhi[0] = 0.0;
		if (pdeg > 0) outdPhi[1] = 1.0;
		for (int i=2;i<=pdeg;++i) {
			outdPhi[i] = TScalar(i)*outPhi[i-1];
		}
		/*for (int i=2;i<=pdeg;++i) {
			outdPhi[i] *= TScalar(i);
		}*/
	}
}



template<typename TExp,typename TScalar>
std::ostream& operator<<(std::ostream& os,const ExponentEngine<TExp,TScalar>& e) {
	os << "Exponent map for p: " << e.GetPdeg() << " base: " << e.GetBaseType() << std::endl;
	for (int i=0;i<e.GetNumShapeFunc();++i) {
		os << "phi[" << i << "] = "
		   << e.GetExpMap()[i].x << ", "
		   << e.GetExpMap()[i].y << ", "
		   << e.GetExpMap()[i].z << std::endl;
	}
	os << "Shape function number: " << e.GetNumShapeFunc() << std::endl;
	os << "end\n";
	return os;
}




#endif /* __Exponent_Engine_hpp__ */
