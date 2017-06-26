#ifndef CELL_HPP_
#define CELL_HPP_

#include "src/models/Element.hpp"

	enum CellType { INNER, CONSTRAINED, BORDER };
	class TriangleCell
	{
	public:
		size_t id;
		size_t nebr[3];
		size_t points[3];
		CellType type;
		point::Point2d c;
		double V;
	public:
		TriangleCell() {};
		TriangleCell(size_t _id) : id(_id) {};
		~TriangleCell() {};
	};

#endif /* CELL_HPP_ */