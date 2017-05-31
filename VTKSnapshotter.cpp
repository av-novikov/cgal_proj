#include "VTKSnapshotter.hpp"

template <class modelType>
const std::string VTKSnapshotter<modelType>::prefix = "";

template<class modelType>
VTKSnapshotter<modelType>::VTKSnapshotter()
{
	pattern = prefix + "GasOil_3D_%{STEP}.vtu";
};
template<class modelType>
VTKSnapshotter<modelType>::~VTKSnapshotter()
{
}
template<class modelType>
std::string VTKSnapshotter<modelType>::getFileName(int i)
{
	std::string filename = pattern;
	return replace(filename, "%{STEP}", std::to_string(i));
};
template<class modelType>
std::string VTKSnapshotter<modelType>::replace(std::string filename, std::string from, std::string to)
{
	size_t start_pos = 0;
	while ((start_pos = filename.find(from, start_pos)) != std::string::npos)
	{
		filename.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return filename;
};
template<class modelType>
void VTKSnapshotter<modelType>::dump(const int i)
{
}