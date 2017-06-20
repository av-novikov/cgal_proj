#ifndef ELEMENT_HPP_
#define ELEMENT_HPP_

#include <initializer_list>
#include <iostream>
#include <array>
#include <utility>

#define EQ_TOL 1.E-8

namespace point
{
	class Point2dContainer
	{
	public:
		static const int dim = 2;
		union
		{
			double coords[dim];
			struct
			{
				double x;		double y;
			};
		};
	};
	class Point3dContainer
	{
	public:
		static const int dim = 3;
		union
		{
			double coords[dim];
			struct
			{
				double x;		double y;		double z;
			};
		};
	};

	template <class Container>
	class AbstractPoint : public Container
	{
	public:
		size_t id;
	public:
		static const int dim = Container::dim;

		AbstractPoint() {};
		AbstractPoint(const AbstractPoint<Container>& m)
		{
			(*this) = m;
		};
		AbstractPoint(std::initializer_list<double> values) : AbstractPoint()
		{
			int i = 0;
			for (auto value : values)
				this->coords[i++] = value;
		};
		AbstractPoint<Container>& operator=(const AbstractPoint<Container>& m)
		{
			for (int i = 0; i < dim; i++)
				this->coords[i] = m.coords[i];

			return *this;
		};
		double norm() const 
		{
			double sum = 0.0;
			for (int i = 0; i < dim; i++)
				sum += this->coords[i] * this->coords[i];
			return sqrt(sum);
		};
	};
	
	template<class Container>
	AbstractPoint<Container> operator-(const AbstractPoint<Container>& m)
	{
		AbstractPoint<Container> result;

		for (int i = 0; i < Container::dim; i++)
			result.coords = -m.coords[i];

		return result;
	}
	template<class Container1, class Container2, class Container3>
	AbstractPoint<Container3> operator+(const AbstractPoint<Container1>& m1, const AbstractPoint<Container2>& m2)
	{
		AbstractPoint<Container3> result;

		for (int i = 0; i < Container3::dim; i++)
			result.coords[i] = m1.coords[i] + m2.coords[i];

		return result;
	}
	template<class Container>
	AbstractPoint<Container> operator+(const AbstractPoint<Container>& m1, const AbstractPoint<Container>& m2)
	{
		return operator+<Container, Container, Container>(m1, m2);
	}
	template<class Container1, class Container2, class Container3>
	AbstractPoint<Container3> operator-(const AbstractPoint<Container1>& m1, const AbstractPoint<Container2>& m2)
	{
		AbstractPoint<Container3> result;

		for (int i = 0; i < Container3::dim; i++)
			result.coords[i] = m1.coords[i] - m2.coords[i];

		return result;
	}
	template<class Container>
	AbstractPoint<Container> operator-(const AbstractPoint<Container>& m1, const AbstractPoint<Container>& m2)
	{
		return operator-<Container, Container, Container>(m1, m2);
	}
	template<class Container>
	AbstractPoint<Container> operator*(const AbstractPoint<Container>& m, const double x)
	{
		AbstractPoint<Container> result;

		for (int i = 0; i < Container::dim; i++)
			result.coords[i] = m.coords[i] * x;

		return result;
	}
	template<class Container>
	AbstractPoint<Container> operator*(const double x, const AbstractPoint<Container>& m)
	{
		return m * x;
	}
	template<class Container>
	AbstractPoint<Container> operator/(const AbstractPoint<Container>& m, const double x)
	{
		return m * (1 / x);
	}
	template<class Container1, class Container2>
	bool operator==(const AbstractPoint<Container1>& m1, const AbstractPoint<Container2>& m2)
	{
		for (int i = 0; i < Container1::dim; i++)
			if (fabs(m1.coords[i] - m2.coords[i]) > EQ_TOL)
				return false;

		return true;
	}
	template<class Container1, class Container2>
	bool operator!=(const AbstractPoint<Container1>& m1, const AbstractPoint<Container2>& m2)
	{
		return !(m1 == m2);
	}
	template<class Container>
	inline double distance(const AbstractPoint<Container>& p1, const AbstractPoint<Container>& p2)
	{
		return (p2 - p1).norm();
	};

	typedef AbstractPoint<Point2dContainer> Point2d;
	typedef AbstractPoint<Point3dContainer> Point3d;
};
namespace std 
{
	template<class Container> 
	inline std::ostream& operator<<(std::ostream& os, const point::AbstractPoint<Container>& pt) 
	{

		os << std::endl;
		for (int i = 0; i < point::AbstractPoint<Container>::dim; i++)
				os << pt.coords[i] << "\t";

		os << std::endl;
		return os;
	}
};
namespace elem
{
	/*enum ELTYPE {EDGE = 2, TRI = 3, TETR = 4, QUAD = 4, HEX = 8};

	template <int N, class TPoint, class TFacet>
	class Element
	{
	public:
		typedef TPoint Point;
		typedef TFacet Facet;
	public:
		static const int size = N;
	public:
		Point c;
		std::array<Point*, N> pts;
		std::array<Facet*, N> fts;

		double V;
	};

	typedef Element<TRI, point::Point2d, Edge<point::Point2d> > Triangle;
	typedef Element<TETR, point::Point3d, Triangle> Tetrahedron;*/

	struct Edge
	{
		double length;
		point::Point2d center;
		size_t nebrInd;
	};
};

#endif /* ELEMENT_HPP_ */