#ifndef SCENE_HPP_
#define SCENE_HPP_

#include <memory>
#include "src/mesh/TriangleMesh.hpp"
#include "paralution.hpp"

template <class modelType, class solverType, class propsType>
class Scene
{
public:
	typedef modelType Model;
	typedef solverType Method;
	typedef propsType Properties;
	typedef typename Model::Mesh Mesh;
protected:
	std::shared_ptr<Model> model;
	std::shared_ptr<Method> method;
public:
	Scene() {};
	~Scene() { paralution::stop_paralution(); };

	void load(const propsType& props, const Task& task)
	{
		model = std::make_shared<Model>();
		model->load(task, props);
		model->setSnapshotter(model.get());

		paralution::init_paralution();

		method = std::make_shared<Method>(model.get());
		model->snapshot_all(0);
	}
};

#endif /* SCENE_HPP_ */