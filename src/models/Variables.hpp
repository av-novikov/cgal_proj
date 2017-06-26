#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

#include <valarray>

namespace var
{
	struct Var1phase
	{
		static const int size = 1;
		typedef std::valarray<double> DataVector;
	};
}

#endif VARIABLES_HPP_