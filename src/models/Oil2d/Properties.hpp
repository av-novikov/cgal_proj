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

		std::vector<double> perm_eff;
		std::vector<double> r_eff;
		std::vector<double> skin;

		double height;
		double p_init;
	};
	struct Properties
	{
		std::vector<double> timePeriods;
		std::vector<double> rates;
		std::vector<double> pwf;
		std::vector<double> skins;
		std::vector<double> radius;

		// If left boundary condition would be 2nd type
		bool leftBoundIsRate;

		std::vector<std::pair<int, int> > perfIntervals;
		std::vector<Skeleton_Props> props_sk;

		double ht;
		double ht_min;
		double ht_max;

		double alpha;

		double r_w;
		double r_e;
		double height;
		double m;
		double kx, ky;
		double dens_sk_stc;
		double beta_sk;

		double visc_oil;
		double dens_oil_stc;
		double beta_oil;
		double b_oil_bore;

		double p_init;
	};
	struct Oil_Props
	{
		double visc;
		double dens_stc;
		double beta;
	};
};

#endif /* OIL2D_PROPERTIES_HPP_ */
