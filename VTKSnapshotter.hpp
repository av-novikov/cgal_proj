#ifndef VTKSNAPSHOTTER_HPP_
#define VTKSNAPSHOTTER_HPP_

#include <string>

template<class modelType>
class VTKSnapshotter
{
public:
	typedef modelType Model;
	typedef typename Model::Mesh Mesh;
private:
	static const std::string prefix;
	std::string pattern;

	const Mesh* mesh;

	std::string replace(std::string filename, std::string from, std::string to);
	std::string getFileName(int i);
public:
	VTKSnapshotter(const Mesh* _mesh);
	~VTKSnapshotter();

	void dump(const int i);
};

#endif /* VTKSNAPSHOTTER_HPP_ */