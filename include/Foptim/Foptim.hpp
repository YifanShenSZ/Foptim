#ifndef Foptim_Foptim_hpp
#define Foptim_Foptim_hpp

#include <Foptim/line-search_1st/steepest_descent.hpp>
#include <Foptim/line-search_1st/steepest_descent_verbose.hpp>
#include <Foptim/line-search_1st/CGDY.hpp>
#include <Foptim/line-search_1st/CGDY_verbose.hpp>
#include <Foptim/line-search_1st/CGPR.hpp>
#include <Foptim/line-search_1st/CGPR_verbose.hpp>

#include <Foptim/line-search_2nd/NewtonRaphson.hpp>
#include <Foptim/line-search_2nd/BFGS.hpp>
#include <Foptim/line-search_2nd/LBFGS.hpp>
#include <Foptim/line-search_2nd/LBFGS_verbose.hpp>

#include <Foptim/least-square/trust_region.hpp>
#include <Foptim/least-square/trust_region_verbose.hpp>
#include <Foptim/least-square/Gauss_BFGS.hpp>

#include <Foptim/constraint/ALagrangian_NewtonRaphson.hpp>
#include <Foptim/constraint/ALagrangian_BFGS.hpp>

#endif