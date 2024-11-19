//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	GridListFactory.h													||
//*	@brief	网格数组工厂, 负责生成网格数组											||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include "GlobalData.h"
#include "GridList.h"

namespace zaran
{
  /// @brief 网格生成器基类
  /// @details 根据控制参数生成网格
  /// @return 返回生成的
  /// @note 只负责生成网格，不负责销毁网格，销毁网格由Field 负责
  class GridFactory
  {
  public:
    /// @brief 生成网格数组
    /// @param grid_list 生成的网格数组
    /// @param grid_num 网格数组大小
    virtual void CreateGrid(GridBase **&grid_list, int& grid_num) {};
    virtual ~GridFactory() {}
  };
} // namespace zaran