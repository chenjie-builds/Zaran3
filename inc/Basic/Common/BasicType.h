#pragma once
#include <vector>
#include <memory>
#include <string>
#include<map>
#include<variant>
#include "Eigen/Dense"
namespace zaran
{
	using std::string;
	using std::map;
	using std::variant;
	template<typename T>
	using Array = std::vector<T>;//array
	template<typename T>
	using Ptr = std::shared_ptr<T>;//pointer

	using DArray = std::vector<double>;//double array
	using IArray = std::vector<int>;//int array
	using SArray = std::vector<string>;//string array
	using Matrix = Eigen::MatrixXd;//double matrix
	using DVector = Eigen::VectorXd;//double vector
	using IVector = Eigen::VectorXi;//int vector
	using DVector2D = Eigen::Vector2d;//double vector2d
	using DVector3D = Eigen::Vector3d;//double vector3d
}
