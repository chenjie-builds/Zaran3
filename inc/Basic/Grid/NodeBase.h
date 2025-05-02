/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file NodeBase.h
 * \brief Node Base class.
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
  class NodeBase
  {
  public:
    NodeBase();
    virtual ~NodeBase();
  public:
    // 返回节点总数
    const count_type& GetCount() const;
  protected:
    void SetCount(count_type count);
  private:
    count_type m_count;
  };
} // namespace zaran