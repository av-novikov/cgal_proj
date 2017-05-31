#ifndef VTKSNAPSHOTTER_HPP_
#define VTKSNAPSHOTTER_HPP_

#include <string>

template<class modelType>
class VTKSnapshotter
{
private:
	static const std::string prefix;
	std::string pattern;

	modelType* model;

	std::string replace(std::string filename, std::string from, std::string to);
	std::string getFileName(int i);
public:
	VTKSnapshotter();
	~VTKSnapshotter();

	void dump(const int i);
};

#endif /* VTKSNAPSHOTTER_HPP_ */