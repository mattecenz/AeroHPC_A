#ifndef INITIALIZATION_HPP
#define INITIALIZATION_HPP

#include "data/SolverData.hpp"

#include "utils/macroUtils.hpp"

void initializeData(Real *data, Real initTime){
	ITERATE_DOMAIN_VELOCITY(i,j,k,false)
  		U(data,i,j,k) = domData.
  	ITERATE_DOMAIN_END()
};

#endif //INITIALIZATION_HPP
