#ifndef VTKSNAPSHOTTER_HPP_
#define VTKSNAPSHOTTER_HPP_

#include <string>
#include <memory>

template<class modelType>
class VTKSnapshotter
{
public:
	typedef modelType Model;
	typedef typename Model::Mesh Mesh;
private:
	static const std::string prefix;
	std::string pattern;

	const Model* model;
	const Mesh* mesh;

	std::string replace(std::string filename, std::string from, std::string to);
	std::string getFileName(int i);
public:
	VTKSnapshotter(const Model* _model);
	~VTKSnapshotter();

	void dump(const int i);
};

#endif /* VTKSNAPSHOTTER_HPP_ */