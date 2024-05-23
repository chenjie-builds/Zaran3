#include "NodeBase.h"
using namespace zaran;
NodeBase::NodeBase()
{
}
NodeBase::~NodeBase()
{
}
const int &NodeBase::GetNodeNum() const
{
    return m_node_num;
}
void NodeBase::SetNodeNum(int node_num)
{
    m_node_num = node_num;
}