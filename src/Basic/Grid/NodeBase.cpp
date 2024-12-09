#include "NodeBase.h"
using namespace zaran;
NodeBase::NodeBase()
{
}
NodeBase::~NodeBase()
{
}
const size_t&NodeBase::GetCount() const
{
    return m_count;
}
void NodeBase::SetCount(size_t count)
{
    m_count = count;
}