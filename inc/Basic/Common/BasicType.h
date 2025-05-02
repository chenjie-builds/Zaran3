/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file BasicType.h
 * \brief Some Basic Type Definitions.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include <vector>
#include <memory>
#include <string>
#include<map>
#include<variant>
#include "Eigen/Dense"
#include "span.hpp"
namespace zaran
{
	using std::string;
	using std::map;
	using std::variant;
	using std::shared_ptr;
	using std::unique_ptr;
	using std::make_shared;
	using std::make_unique;
	using nonstd::span_lite::span;

	template<typename T>
	using dynamic_array = std::vector<T>;//array
	//index type
	using index_type = size_t;
	//count type
	using count_type = size_t;//count type

	//dimension type
	using dimension_type = size_t;



}
