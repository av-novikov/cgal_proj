#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

#include <valarray>
#include "adolc/adouble.h"

namespace var
{
	namespace containers
	{
		struct Var1phase
		{
			static const int size = 1;
			double& p;
		};
		struct TapeVar1Phase
		{
			static const int size = 1;
			adouble p;
		};
		struct AcidVar
		{
			static const int size = 5;
			double& m;
			double& p;
			double& s;
			double& xa;
			double& xw;
		};
		struct TapeAcidVar
		{
			static const int size = 5;
			adouble m;
			adouble p;
			adouble s;
			adouble xa;
			adouble xw;
		};
	};

	template <typename TVariable>
	struct BaseVarWrapper
	{
		TVariable u_prev, u_iter, u_next;
	};
	template <typename TVariable>
	struct BasicVariables
	{
		static const int size = TVariable::size;
		typedef BaseVarWrapper<TVariable> Wrap;
		std::valarray<double> u_prev, u_iter, u_next;

		Wrap operator[](const size_t idx)
		{
			return{ { u_prev[idx * size] },{ u_iter[idx * size] },{ u_next[idx * size] } };
		};
		const Wrap operator[](const size_t idx) const
		{
			return{ { u_prev[idx * size] },{ u_iter[idx * size] },{ u_next[idx * size] } };
		};
	};
}

#endif /* VARIABLES_HPP_ */