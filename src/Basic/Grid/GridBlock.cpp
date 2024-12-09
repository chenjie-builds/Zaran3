#include "GridBlock.h"
namespace zaran
{
    GridBlock::GridBlock(const string &name, Id index, Id dim) : GridStruct(name, index, dim)
    {
    }
    GridBlock::~GridBlock()
    {
    }

    void GridBlock::Allocate(Id ni, Id nj, Id nk, Id ghost_level)
    {
        GridStruct::Allocate(ni, nj, nk, ghost_level);
        m_iblank.resize((ni + 2 * ghost_level) * (nj + 2 * ghost_level) * (nk + 2 * ghost_level));
    }
    void GridBlock::SetIBlank(Id i, Id j, Id k, IBlank iblank)
    {
        m_iblank[i + j * GetNi() + k * GetNi() * GetNj()] = iblank;
    }
    const IBlank &GridBlock::GetIBlank(Id i, Id j, Id k) const
    {
        return m_iblank[i + j * GetNi() + k * GetNi() * GetNj()];
    }
}