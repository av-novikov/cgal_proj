#ifndef CELL_HPP_
#define CELL_HPP_

#include "Element.hpp"

namespace cell
{
	template <class TElement, typename TVariable>
	class AbstractCell : public TElement
	{
	public:
		typedef TElement Element;
		typedef typename Element::Point Point;
		typedef typename Element::Facet Facet;
		typedef TVariable Variable;
	public:
		Variable u_prev, u_iter, u_next;
		std::array<AbstractCell*, Element::size> nebrs;
	};

	template <typename TVariable>
	class TriangleCell : public AbstractCell<elem::Triangle, TVariable>
	{
	public:
	};
};

#endif /* CELL_HPP_ */