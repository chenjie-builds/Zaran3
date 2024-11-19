#include "GridBlock.h"
namespace zaran
{
    GridBlock::GridBlock(const string &name, int index, int dim) : GridStruct(name, index, dim)
    {
    }
    GridBlock::~GridBlock()
    {
    }

    void GridBlock::Allocate(int ni, int nj, int nk, int ghost_level)
    {
        GridStruct::Allocate(ni, nj, nk, ghost_level);
        m_iblank = new IBlank[(ni + 2 * ghost_level) * (nj + 2 * ghost_level) * (nk + 2 * ghost_level)];
    }
    void GridBlock::SetIBlank(int i, int j, int k, IBlank iblank)
    {
        m_iblank[i + j * GetNi() + k * GetNi() * GetNj()] = iblank;
    }
    const IBlank &GridBlock::GetIBlank(int i, int j, int k) const
    {
        return m_iblank[i + j * GetNi() + k * GetNi() * GetNj()];
    }
}