#ifndef ACID2D_PROPERTIES_HPP_
#define ACID2D_PROPERTIES_HPP_

#include <array>

#include "src/models/Acid/Reactions.hpp"

namespace acid2d
{
	struct Skeleton_Props : public basic2d::Skeleton_Props
	{
		//SolidComponent cur_mineral;
		double xa_eqbm;

		// Initial values
		double p_init;
		double s_init;
		double xa_init;
		double xw_init;

		//double d_pore_r, d_pore_z;
		inline adouble getPermCoseni_x(adouble m) const
		{
			//return d_pore_r * d_pore_r * m * m * m / (1 - m) / (1 - m) / 150.0;
			return kx * (m * m * m / (1 - m) / (1 - m)) / 
							(m_init * m_init * m_init / (1 - m_init) / (1 - m_init));
		};
		inline adouble getPermCoseni_y(adouble m) const
		{
			//return d_pore_z * d_pore_z * m * m * m / (1 - m) / (1 - m) / 150.0;
			return ky * (m * m * m / (1 - m) / (1 - m)) /
				(m_init * m_init * m_init / (1 - m_init) / (1 - m_init));
		};
		inline double getInitDiam(double m_init, double k0)
		{
			return sqrt(150.0 * k0 * (1.0 - m_init) * (1.0 - m_init) / m_init / m_init / m_init);
		};

		inline adouble getDensity(adouble p) const
		{
			return dens_stc;
		};
	};
	struct Water_Props : public basic2d::Liquid_Props
	{
		//LiquidComponent acid;
		//SolidComponent salt;
		//LiquidComponent water;

		Interpolate* kr;
		inline adouble getKr(adouble s, const Skeleton_Props& props) const
		{
			adouble isAboveZero = (s - props.s_wc > 0.0) ? true : false;
			adouble isAboveCritical = (s > 1.0 - props.s_oc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, pow((s - (adouble)props.s_wc) / (adouble)(1.0 - props.s_wc - props.s_oc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)1.0);
			return tmp;
		};
		inline adouble getViscosity(adouble p, adouble xa, adouble xw) const
		{
			return visc;
		};
		inline adouble getDensity(adouble p, adouble xa, adouble xw) const
		{
			return dens_stc;
		};
	};
	struct Oil_Props : public basic2d::Liquid_Props
	{
		double gas_dens_stc;
		//LiquidComponent oil;
		Interpolate* b;
		inline adouble getB(adouble p) const
		{
			return exp(-(adouble)beta * (p - p_ref));
		};
		inline adouble getDensity(adouble p) const
		{
			return dens_stc / getB(p);
		};

		Interpolate* kr;
		inline adouble getKr(adouble s, const Skeleton_Props& props) const
		{
			adouble isAboveZero = (1.0 - s - props.s_oc > 0.0) ? true : false;
			adouble isAboveCritical = (1.0 - s > 1.0 - props.s_wc) ? true : false;
			adouble tmp;
			condassign(tmp, isAboveZero, pow((1.0 - s - (adouble)props.s_oc) / (adouble)(1.0 - props.s_wc - props.s_oc), 3.0), (adouble)0.0);
			condassign(tmp, isAboveCritical, (adouble)1.0);
			return tmp;
		};
		inline adouble getViscosity(adouble p) const
		{
			return visc;
		};
	};
	struct Properties : public basic2d::Properties
	{
		std::vector<Skeleton_Props> props_sk;
		Water_Props props_w;
		Oil_Props props_o;

		std::vector<double> xa;

		std::vector< std::pair<double, double> > rho_co2;
	};
};

#endif /* ACID2D_PROPERTIES_HPP_ */