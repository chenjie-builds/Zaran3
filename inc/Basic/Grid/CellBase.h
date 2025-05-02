/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file CellBase.h
 * \brief Cell Topology Base class.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
namespace zaran
{
class CellBase
{
  public:
    CellBase();
    virtual ~CellBase();
    const int &GetCellNum() const;

  protected:
    void SetCellNum(int cell_num);
  private:
    int m_cell_num;
};

} // namespace zaran