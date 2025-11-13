/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file CommonPara.h
 * \brief Define some	 constant parameter for global use
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
namespace zaran
{
// 空间维数
enum class Dimension
{
    one = 1,
    two = 2,
    three = 3
};
const dimension_type ONE_DIM = 1;
const dimension_type TWO_DIM = 2;
const dimension_type THREE_DIM = 3;

const index_type ID_DENSITY = 0;
const index_type ID_VELOCITY_X = 1;
const index_type ID_VELOCITY_Y = 2;
const index_type ID_VELOCITY_Z = 3;
const index_type ID_PRESSURE = 4;
const index_type ID_TEMPERATURE = 5;

const double LARGE_NUMBER = 1.0e40;
const double SMALL_NUMBER = 1.0e-40;
const double TINY_NUMBER = 1.0e-40;
const double EPSILON_NUMBER = 1.0e-12;

const double PI = 3.141592653589793238462643383279502884197169399; // 圆周率pi
struct Box
{
    double x_min = LARGE_NUMBER;
    double x_max = -LARGE_NUMBER;
    double y_min = LARGE_NUMBER;
    double y_max = -LARGE_NUMBER;
    double z_min = LARGE_NUMBER;
    double z_max = -LARGE_NUMBER;
};
} // namespace zaran