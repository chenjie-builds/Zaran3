#include "GridFNFDM.h"
namespace zaran
{
    GridFN::GridFN(const string& name, int index, dimension_type dim) :GridBase(name, index, dim, GridType::Flexible)
    {
		m_node = make_unique<NodeFN>(0);
        m_face = make_unique<FaceFN>(0);
		m_cell = make_unique<CellFN>(0);
		m_boundary_map = make_unique<BoundManagerFN>();
		m_inner_node.resize(0);
		m_bound_node_index.resize(0);
    }
    GridFN::~GridFN()
    {

    }
	int GridFN::GetTotalNodeNum() const
    {
		return m_node->GetCount();
    }
    int GridFN::GetInnerNodeNum() const
    {
		return m_inner_node.size();
    }
    int GridFN::GetBoundNodeNum() const
    {
        return m_bound_node_index.size();
    }

	NodeFN& GridFN::GetNode() 
	{
        return *m_node;
	}

	FaceFN& GridFN::GetFace()
    {
        return *m_face;
    }
    CellFN& GridFN::GetCell()
    {
        return *m_cell;
    }
    BoundManagerFN& GridFN::GetBoundaryMap()
    {
        return *m_boundary_map;
    }
    void GridFN::InitNode()
    {
        count_type node_count = m_node->GetCount();
		count_type inner_node_count = 0;
		count_type bound_node_count = 0;
        for (int iNode = 0; iNode < node_count; iNode++)
        {
            if (m_node->GetType(iNode) == NodeType::inner)
            {
                inner_node_count++;
            }
            else
            {
                bound_node_count++;
            }
        }
		m_inner_node.resize(inner_node_count);
		m_bound_node_index.resize(bound_node_count);
        int inner_index = 0;
        int bound_index = 0;
        for (int iNode = 0; iNode < node_count; iNode++)
        {
            if (m_node->GetType(iNode) == NodeType::inner)
            {
                m_inner_node[inner_index++] = iNode;
            }
            else
            {
                m_bound_node_index[bound_index++] = iNode;
            }
        }
    }
}
