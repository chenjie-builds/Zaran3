#include "IdxProxy.h"
namespace zaran
{
    Idx::Idx(int idx)
    {
        m_idx = idx;
    }

    Idx::~Idx()
    {
    }
    int Idx::GetIdx() const
    {
        return m_idx;
    }
    void Idx::SetIdx(int idx)
    {
        m_idx = idx;
    }
}