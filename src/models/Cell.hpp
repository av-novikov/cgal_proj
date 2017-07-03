#ifndef CELL_HPP_
#define CELL_HPP_

#include "src/models/Element.hpp"

	enum CellType { INNER, FRAC, BORDER, WELL };
	class TriangleCell
	{
	public:
		size_t id;
		size_t nebr[3];
		double dist[3];
		double length[3];
		size_t points[3];
		CellType type;
		point::Point2d c;
		double V;
	public:
		TriangleCell() {};
		TriangleCell(size_t _id) : id(_id) {};
		~TriangleCell() {};

		const double& getDistance(const size_t idx) const 
		{
			for (int i = 0; i < 3; i++)
				if (nebr[i] == idx)
					return dist[i];
		};
	};

#endif /* CELL_HPP_ */