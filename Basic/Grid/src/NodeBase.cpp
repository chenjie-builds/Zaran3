#include "NodeBase.h"
using namespace zaran;
NodeBase::NodeBase()
{
}
NodeBase::~NodeBase()
{
}
const int &NodeBase::GetCount() const
{
    return m_count;
}
void NodeBase::SetCount(int count)
{
    m_count = count;
}