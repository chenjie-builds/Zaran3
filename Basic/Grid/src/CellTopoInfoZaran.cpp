#include "CellTopoInfoZaran.h"
namespace zaran
{

    CellTopoInfoZaran::CellTopoInfoZaran()
    {
        CellTopoInfo::CellTopoInfo();
        m_cellType.resize(0);
    }

    CellTopoInfoZaran::~CellTopoInfoZaran()
    {
    }

    void CellTopoInfoZaran::SetType(Array<CellType>& cellType)
    {
        m_cellType = cellType;
    }

    Array<CellType>& CellTopoInfoZaran::GetType()
    {
        return m_cellType;
    }
}