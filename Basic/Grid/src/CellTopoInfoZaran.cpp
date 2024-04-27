#include "CellTopoInfoZaran.h"
namespace zaran
{

    CellTopoZaran::CellTopoZaran()
    {
        CellTopo::CellTopo();
        m_cellType.resize(0);
    }

    CellTopoZaran::~CellTopoZaran()
    {
    }

    void CellTopoZaran::SetType(Array<CellType>& cellType)
    {
        m_cellType = cellType;
    }

    Array<CellType>& CellTopoZaran::GetType()
    {
        return m_cellType;
    }
}