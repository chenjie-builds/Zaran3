/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file CellStruct.h
 * \brief Cell Topology Class for Structured Grid.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "CellBase.h"
namespace zaran
{
class CellStruct : public CellBase
{
  public:
    CellStruct() ;
    virtual ~CellStruct();
    void Allocate(int i_num, int j_num, int k_num);
    void SetCenterCoord(int idx_i, int idx_j, int idx_k, double *center);
    const double *GetCenterCoord(int idx_i, int idx_j, int idx_k) const;
  protected:
    int GetIdx(int idx_i, int idx_j, int idx_k) const;
  private:
    /// @brief 单元个数，分别为i、j、k方向的单元个数，显然每个方向的单元个数是节点个数减一
    int m_i_num, m_j_num, m_k_num;
    /// @brief 中心坐标
    double *m_center_coord;
};
} // namespace zaran