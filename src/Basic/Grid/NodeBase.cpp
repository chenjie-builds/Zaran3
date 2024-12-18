#include "NodeBase.h"
using namespace zaran;
NodeBase::NodeBase()
{
}
NodeBase::~NodeBase()
{
}
const count_type&NodeBase::GetCount() const
{
    return m_count;
}
void NodeBase::SetCount(count_type count)
{
    m_count = count;
}