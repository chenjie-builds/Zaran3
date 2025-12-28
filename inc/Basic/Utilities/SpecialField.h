/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file SpecialField.h
 * \brief Utility functions for generating special flow fields.
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
void CalcSupersonicVortex(const double &x, const double &y, const double &x0, const double y0, const double &ci,
                          const double &Mi, const double &rhoi, Eigen::VectorXd &primitive);
void CalcIsentropicVortex(const double &x, const double &y, const double &beta, Eigen::VectorXd &primitive);
} // namespace zaran