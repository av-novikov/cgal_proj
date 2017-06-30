#ifndef OIL2D_PROPERTIES_HPP_
#define OIL2D_PROPERTIES_HPP_

#include <vector>
#include <utility>

namespace oil2d
{
	struct Skeleton_Props
	{
		double m;
		double beta;
		double dens_stc;

		double kx, ky;

		double height;
		double p_init;
		double p_out;
	};
	struct Oil_Props
	{
		double visc;
		double dens_stc;
		double beta;
	};
	struct Properties
	{
		std::vector<double> timePeriods;
		std::vector<double> rates;
		std::vector<double> pwf;

		bool leftBoundIsRate;
		bool rightBoundIsPres;

		std::vector<std::pair<int, int> > perfIntervals;
		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;

		double ht;
		double ht_min;
		double ht_max;

		double alpha;

		double r_w;
		double r_e;
	};
};

#endif /* OIL2D_PROPERTIES_HPP_ */
