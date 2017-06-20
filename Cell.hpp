#ifndef CELL_HPP_
#define CELL_HPP_

#include "Element.hpp"

	enum CellType { INNER, BORDER1, BORDER2, BORDER3, CONSTRAINED };
	class TriangleCell
	{
	public:
		size_t id;
		size_t nebr[3];
		CellType type;
		point::Point2d c;
		double V;
	public:
		TriangleCell() {};
		TriangleCell(size_t _id) : id(_id) {};
		~TriangleCell() {};
	};

#endif /* CELL_HPP_ */