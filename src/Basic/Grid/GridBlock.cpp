#include "GridBlock.h"
namespace zaran
{
    GridBlock::GridBlock(const string &name, index_type index, index_type dim) : GridStruct(name, index, dim)
    {
    }
    GridBlock::~GridBlock()
    {
    }

    void GridBlock::Allocate(index_type ni, index_type nj, index_type nk, index_type ghost_level)
    {
        GridStruct::Allocate(ni, nj, nk, ghost_level);
        m_iblank.resize((ni + 2 * ghost_level) * (nj + 2 * ghost_level) * (nk + 2 * ghost_level));
    }
    void GridBlock::SetIBlank(index_type i, index_type j, index_type k, IBlank iblank)
    {
        m_iblank[i + j * GetNi() + k * GetNi() * GetNj()] = iblank;
    }
    const IBlank &GridBlock::GetIBlank(index_type i, index_type j, index_type k) const
    {
        return m_iblank[i + j * GetNi() + k * GetNi() * GetNj()];
    }
}