#ifndef INITIALIZATION_HPP
#define INITIALIZATION_HPP

#include "data/SolverData.hpp"

#include "utils/macroUtils.hpp"

inline void initializeModel(Real *data, const Real initTime){
	{ITERATE_DOMAIN_VELOCITY(i,j,k,true, NO_SKIP)
  		U(data,i,j,k) = domData.domainF.UF(x,y,z,initTime);
  		V(data,i,j,k) = domData.domainF.VF(x,y,z,initTime);
  		W(data,i,j,k) = domData.domainF.WF(x,y,z,initTime);
  	ITERATE_DOMAIN_END()}
    {ITERATE_DOMAIN_PRESSURE(i,j,k,true)
    	P(data,i,j,k) = domData.domainF.PF(x,y,z,initTime);
    ITERATE_DOMAIN_END()}
};

#endif //INITIALIZATION_HPP
