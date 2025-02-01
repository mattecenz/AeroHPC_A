#ifndef INITIALIZATION_HPP
#define INITIALIZATION_HPP

#include "data/SolverData.hpp"

#include "utils/macroUtils.hpp"

inline void initializeModel(Real *data, const Real initTime){
	ITERATE_DOMAIN_VELOCITY(i,j,k,true, NO_SKIP)
  		U(data,i,j,k) = domData.domainF.UF(X,Y,Z,initTime);
  		V(data,i,j,k) = domData.domainF.VF(X,Y,Z,initTime);
  		W(data,i,j,k) = domData.domainF.WF(X,Y,Z,initTime);
  	ITERATE_DOMAIN_END()
    ITERATE_DOMAIN_PRESSURE(i,j,k,true)
    	P(data,i,j,k) = domData.domainF.PF(X,Y,Z,initTime);
    ITERATE_DOMAIN_END()
};

#endif //INITIALIZATION_HPP
