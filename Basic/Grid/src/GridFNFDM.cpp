#include "GridFNFDM.h"
namespace zaran
{
    GridFN::GridFN(const string& name, int index, int dim) :GridBase(name, index, dim, GridType::Flexible)
    {
        m_node = nullptr;
        m_face = nullptr;
        m_cell = nullptr;
        m_boundary_map = nullptr;
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
    }
    void GridFN::SetNode(NodeFN* node)
    {
        m_node = node;
    }
    void GridFN::SetFace(FaceFN* face)
    {
        m_face = face;
    }
    void GridFN::SetCell(CellFN* cell)
    {
        m_cell = cell;
    }
    void GridFN::SetBoundaryMap(BoundaryMap* boundaryMap)
    {
        m_boundary_map = boundaryMap;
    }
    int GridFN::GetTotalNodeNum() const
    {
        return m_node->GetNodeNum();
    }
    int GridFN::GetInnerNodeNum() const
    {
        return 0;
    }
    int GridFN::GetBoundNodeNum() const
    {
        return 0;
    }
    NodeFN* GridFN::GetNodeTopo()
    {
        return m_node;
    }
    FaceFN* GridFN::GetFaceTopo()
    {
        return m_face;
    }
    CellFN* GridFN::GetCellTopo()
    {
        return m_cell;
    }
    BoundaryMap* GridFN::GetBoundaryMap()
    {
        return m_boundary_map;
    }
} // namespace zaran