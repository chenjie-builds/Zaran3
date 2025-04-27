#include "GridStruct.h"
namespace zaran
{
	GridStruct::GridStruct(const string &name, index_type index, index_type dim) : GridBase(name, index, dim, GridType::Structured)
	{
		m_ghost_level = 1;
	}
	GridStruct::~GridStruct()
	{
	}
	void GridStruct::Allocate(index_type ni, index_type nj, index_type nk, index_type ghost_level)
	{
		m_ghost_level = ghost_level;
		m_node = make_unique<NodeStruct>(ni + ghost_level * 2, nj + ghost_level * 2, nk + ghost_level * 2);
		m_face = make_unique<FaceStruct>();
		m_cell = make_unique<CellStruct>();
		m_cell->Allocate(ni + ghost_level * 2 - 1, nj + ghost_level * 2 - 1, nk + ghost_level * 2 - 1);
		m_bound_map = make_unique<BoundStructManager>();
		m_idx_proxy = make_shared<IdProxyStruct>(ni + ghost_level * 2, nj + ghost_level * 2, nk + ghost_level * 2);
	}
	int GridStruct::GetNi()const
	{
		return m_node->GetINum();
	}

	int GridStruct::GetNj()const
	{
		return m_node->GetJNum();
	}

	int GridStruct::GetNk()const
	{
		return m_node->GetKNum();
	}
	void GridStruct::GetNodeNum(index_type&ni, index_type&nj, index_type&nk) const
	{
		ni = GetNi();
		nj = GetNj();
		nk = GetNk();
	}
	int GridStruct::GetTotalNodeNum()
	{
		return m_node->GetCount();
	}
	void GridStruct::GetRange(int &is, int &ie, int &js, int &je, int &ks, int &ke) const
	{
		is = GetGhostLevel();
		ie = GetNi() - GetGhostLevel() - 1;
		js = GetGhostLevel();
		je = GetNj() - GetGhostLevel() - 1;
		ks = GetGhostLevel();
		ke = GetNk() - GetGhostLevel() - 1;
	}
}