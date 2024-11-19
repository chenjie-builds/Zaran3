#include "IdxProxy.h"
namespace zaran
{
    IdxProxy::IdxProxy(int idx)
    {
        m_idx = idx;
    }

    IdxProxy::~IdxProxy()
    {
    }
    int IdxProxy::GetIdx() const
    {
        return m_idx;
    }
    void IdxProxy::SetIdx(int idx)
    {
        m_idx = idx;
    }
}