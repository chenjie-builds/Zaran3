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
    virtual GridBase* CreateGrid() = 0;
    virtual ~GridFactory() {}
};
} // namespace zaran