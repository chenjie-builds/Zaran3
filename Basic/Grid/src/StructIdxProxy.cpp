#include "StructIdxProxy.h"
namespace zaran
{
    StructIdxProxy::StructIdxProxy(GridStruct* grid) :m_grid(grid)
    {
    }
    StructIdxProxy::~StructIdxProxy()
    {

    }
    int StructIdxProxy::GetNi() const
    {
        return m_grid->GetNi();
    }
    int StructIdxProxy::GetNj() const
    {
        return m_grid->GetNj();
    }
    int StructIdxProxy::GetNk() const
    {
        return m_grid->GetNk();
    }
    int StructIdxProxy::GetIdx(int i, int j, int k)const
    {
        return i + j * m_grid->GetNi() + k * m_grid->GetNi() * m_grid->GetNj();
    }
    void StructIdxProxy::GetIdxStruct(int idx, int& i, int& j, int& k) const
    {
        int ni = m_grid->GetNi();
        int nj = m_grid->GetNj();
        i = idx % ni;
        j = (idx / ni) % nj;
        k = idx / (ni * nj);
    }
}