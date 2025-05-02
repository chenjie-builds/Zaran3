/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GridGenerator.h
 * \brief Grid factory, Generate grids.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */

#pragma once
#include "GlobalData.h"
#include "GridBase.h"
namespace zaran
{
  /// @brief 网格生成器基类
  /// @details 根据控制参数生成网格
  /// @return 返回生成的
  /// @note 只负责生成网格，不负责销毁网格，销毁网格由Field 负责
  class GridGenerator
  {
  public:
    /// @brief 生成网格数组
    /// @param grid_list 生成的网格数组
    /// @param grid_num 网格数组大小
    virtual void CreateGrid(dynamic_array<shared_ptr<GridBase>>&grid_list) {};
    virtual ~GridGenerator() {}
  };
} // namespace zaran