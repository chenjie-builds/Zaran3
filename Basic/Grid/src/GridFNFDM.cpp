#include "GridFNFDM.h"
namespace zaran
{
    GridFN::GridFN(const string& name, int index, int dim) :GridBase(name, index, dim, GridType::Flexible)
    {
        m_node = nullptr;
        m_face = nullptr;
        m_cell = nullptr;
        m_boundary_map = nullptr;
        m_inner_node_num = 0;
        m_inner_node_index = nullptr;
        m_bound_node_num = 0;
        m_bound_node_index = nullptr;
        m_total_node_num = 0;
    }
    GridFN::~GridFN()
    {
        if (m_node)
        {
            delete m_node;
            m_node = nullptr;
        }
        if (m_face)
        {
            delete m_face;
            m_face = nullptr;
        }
        if (m_cell)
        {
            delete m_cell;
            m_cell = nullptr;
        }
        if (m_boundary_map)
        {
            delete m_boundary_map;
            m_boundary_map = nullptr;
        }
        if (m_inner_node_index)
        {
            delete[] m_inner_node_index;
            m_inner_node_index = nullptr;
        }
        if (m_bound_node_index)
        {
            delete[] m_bound_node_index;
            m_bound_node_index = nullptr;
        }
    }
    void GridFN::SetNode(NodeFN* node)
    {
        m_node = node;
        if (m_node)
            InitNode();
    }
    void GridFN::SetFace(FaceFN* face)
    {
        m_face = face;
    }
    void GridFN::SetCell(CellFN* cell)
    {
        m_cell = cell;
    }
    void GridFN::SetBoundaryMap(BoundMapFN* boundaryMap)
    {
        m_boundary_map = boundaryMap;
    }
    int GridFN::GetTotalNodeNum() const
    {
        return m_total_node_num;
    }
    int GridFN::GetInnerNodeNum() const
    {
        return m_inner_node_num;
    }
    int GridFN::GetBoundNodeNum() const
    {
        return m_bound_node_num;
    }
    NodeFN* GridFN::GetNode()
    {
        return m_node;
    }
    FaceFN* GridFN::GetFace()
    {
        return m_face;
    }
    CellFN* GridFN::GetCell()
    {
        return m_cell;
    }
    BoundMapFN* GridFN::GetBoundaryMap()
    {
        return m_boundary_map;
    }
    int* GridFN::GetInnerNode()
    {
        return m_inner_node_index;
    }
    int* GridFN::GetBoundNode()
    {
        return m_bound_node_index;
    }
    void GridFN::InitNode()
    {
        m_total_node_num = m_node->GetNodeNum();
        for (int iNode = 0; iNode < m_total_node_num; iNode++)
        {
            if (m_node->GetType(iNode) == NodeType::inner)
            {
                m_inner_node_num++;
            }
            else
            {
                m_bound_node_num++;
            }
        }
        m_inner_node_index = new int[m_inner_node_num];
        m_bound_node_index = new int[m_bound_node_num];
        int inner_index = 0;
        int bound_index = 0;
        for (int iNode = 0; iNode < m_total_node_num; iNode++)
        {
            if (m_node->GetType(iNode) == NodeType::inner)
            {
                m_inner_node_index[inner_index] = iNode;
                inner_index++;
            }
            else
            {
                m_bound_node_index[bound_index] = iNode;
                bound_index++;
            }
        }
    }
}
