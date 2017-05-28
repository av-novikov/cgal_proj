#ifndef VARIABLES_HPP_
#define VARIABLES_HPP_

namespace var
{
	struct Var1phase
	{
		static const int size = 1;
		union {
			double values[1];
			struct
			{
				double p;
			};
		};
	};
}

#endif VARIABLES_HPP_