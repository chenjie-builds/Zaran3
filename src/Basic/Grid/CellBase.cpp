#include "CellBase.h"
namespace zaran
{
CellBase::CellBase()
{
}
CellBase::~CellBase()
{
}
const int &CellBase::GetCellNum() const
{
    return m_cell_num;
}
void CellBase::SetCellNum(int cell_num)
{
    m_cell_num = cell_num;
}
} // namespace zaran