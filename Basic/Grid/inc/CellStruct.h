#pragma once
#include "CellBase.h"
namespace zaran
{
class CellStruct : public CellBase
{
  public:
    CellStruct(int ni, int nj, int nk) ;
    virtual ~CellStruct();
    void SetCenterCoord(int i, int j, int k, double *center);
    const double *GetCenterCoord(int i, int j, int k) const;
  protected:
    int GetIndex(int i, int j, int k) const;
  private:
    int m_ni, m_nj, m_nk;
    double *m_center;
};
} // namespace zaran