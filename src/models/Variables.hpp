#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

#include <valarray>

namespace var
{
	namespace containers
	{
		struct Var1phase
		{
			static const int size = 1;
			double& p;
		};
		struct Vat2phase
		{
			static const int size = 2;
			double& p;
			double& s;
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