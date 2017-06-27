#ifndef ABSTRACTMODEL_HPP_
#define ABSTRACTMODEL_HPP_

#include <vector>
#include <string>
#include <map>
#include <memory>

#include "src/snapshotter/VTKSnapshotter.hpp"
#include "src/mesh/TriangleMesh.hpp"
//#include "util/utils.h"

template <typename TVarContainer, typename propsType, template <typename TVarContainer> class TVariables>
class AbstractModel : public TVariables
{	
public:
	typedef TVarContainer VarContainer;
	typedef TVariables Variables;
	typedef mesh::TriangleMesh<VarContainer> Mesh;
protected:
	std::shared_ptr<Mesh> mesh;

	// Spacial properties
	double length_perf;
	double r_w;
	double r_e;
	double Volume;
	int cellsNum;
		
	// Rate of the well
	double Q_sum;
	double Pwf;
	// Ranges of perforated cells numbers
	std::vector<std::pair<int,int> > perfIntervals;
	// Vector of <cell number, rate in the cell> for left border cells
	std::map<int,double> Qcell;

	// Temporary properties
	double ht;
	double ht_min;
	double ht_max;
		
	// Number of periods
	int periodsNum;
	// End times of periods [sec]
	std::vector<double> period;
	// Oil rates [m3/day]
	std::vector<double> rate;
	// Vector of BHPs [bar]
	std::vector<double> pwf;
	// If left boundary condition would be 2nd type
	bool leftBoundIsRate;
	// If right boundary condition would be 1st type
	bool rightBoundIsPres;
	// BHP will be converted to the depth
	double depth_point;
	// During the time flow rate decreases 'e' times in well test [sec] 
	double alpha;
	double wellboreDuration;

	virtual void loadMesh() = 0;
	virtual void setProps(propsType& props) = 0;
	virtual void makeDimLess() = 0;
	virtual void setPerforated() = 0;
	virtual void setInitialState() = 0;

	static const int var_size;
public:
	AbstractModel() { isWriteSnaps = true;	grav = 9.8; };
	virtual ~AbstractModel() {};
	
	// Dimensions
	double t_dim;
	double R_dim;
	double Z_dim;
	double P_dim;
	double T_dim;
	double Q_dim;
	double grav;

	void load(propsType& props)
	{
		setProps(props);

		loadMesh();
		setPerforated();
		setInitialState();
	};
	virtual void setPeriod(int period) = 0;
	virtual void setWellborePeriod(int period, double cur_t) {};
	int getCellsNum() {	return cellsNum; };

	void snapshot_all(int i) { snapshotter->dump_all(i); }
};

#endif /* ABSTRACTMODEL_HPP_ */
