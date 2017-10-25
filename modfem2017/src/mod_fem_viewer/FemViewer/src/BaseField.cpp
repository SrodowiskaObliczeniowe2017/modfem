/*
 * BaseField.cpp
 *
 *  Created on: 2011-06-03
 *      Author: pawel
 */
#include "Log.h"
#include "fv_config.h"
#include "aph_intf.h"
#include "BaseField.h"

#ifdef _USE_FV_EXT_MOD
#include "ApproxModule.h"
#endif

namespace FemViewer
{

const char FieldType[] = {'s','c','t','m'};

#ifdef _USE_FV_LIB
int_global_solutiom BaseField::apr_global_solution = apr_sol_xglob;
double_el_calc BaseField::apr_elem_calc_3D_mod = apr_elem_calc_3D;
intfcharp BaseField::module_introduce = apr_module_introduce;
intfchar5intcharpfnp BaseField::init_field = apr_init_field;
intfint BaseField::free_field = apr_free_field;
intf2int BaseField::get_base_type = apr_get_base_type;
intfint BaseField::get_nr_sol = apr_get_nr_sol;
intfint BaseField::get_nreq = apr_get_nreq;
intf2intintp BaseField::get_el_pdeg = apr_get_el_pdeg;
intf2intintp BaseField::get_el_pdeg_numshap = apr_get_el_pdeg_numshap;
intf3intdoublep BaseField::get_element_dofs = apr_get_el_dofs;
#else
int_global_solutiom BaseField::apr_global_solution = NULL;
double_el_calc BaseField::apr_elem_calc_3D_mod = NULL;
intfcharp BaseField::module_introduce = NULL;
intfchar5intcharpfnp BaseField::init_field = NULL;
intfint BaseField::free_field = NULL;
intf2int BaseField::get_base_type = NULL;
intfint BaseField::get_nr_sol = NULL;
intfint BaseField::get_nreq = NULL;
intf2intintp BaseField::get_el_pdeg = NULL;
intf2intintp BaseField::get_el_pdeg_numshap = NULL;
intf3intdoublep BaseField::get_element_dofs = NULL;
#endif

int BaseField::GetApproximationType()
{
	//mfp_debug("BaseField::GetApproximationType\n");
	char type[32];
	module_introduce(type);
	int result;
	if (type[0] == 'D') result = FieldDG;
	else if(type[0] =='S') result = FieldSTD;
	else result = Unknown;
	//mfp_debug("After introduce %d\n",result);
	return result;
}


int BaseField::Init(const char* fname_,const int parent_id_)
{

	#ifdef _USE_FV_LIB
	Id()  = APC_CUR_FIELD_ID;
	#else
    char type    = FieldType[this->_htype];
    Id()  = init_field(type,APC_READ,parent_id_,-1,-1,-1,const_cast<char*>(fname_),NULL);
    if (Id() < 0) return(-1);
	#endif
	GetType() = static_cast<HandleType>(this->GetApproximationType());
    _solution.nrSolutions = get_nr_sol(this->idx());
    _solution.nrEquetions = get_nreq(this->idx());

	return(Id());
}

int BaseField::Free()
{
	#ifndef _USE_FV_LIB
	return free_field(this->idx());
	#else
	return 0;
	#endif
}

//int BaseField::Reload()
//{
// return 0;
//}
//
//int BaseField::Update()
//{
//  return 0;
//}

void BaseField::Reset()
{
	Free();
	BaseHandle::Reset();
	_solution.Clear();
}

}// end namespace FemViewer
