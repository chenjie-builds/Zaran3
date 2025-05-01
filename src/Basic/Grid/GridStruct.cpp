#include "GridStruct.h"
namespace zaran
{
	GridStruct::GridStruct(const string& name, index_type index, index_type dim) : GridBase(name, index, dim, GridType::Structured)
	{
		m_ghost_level = 1;
	}
	GridStruct::~GridStruct()
	{
	}
	void GridStruct::Allocate(count_type ni, count_type nj, count_type nk, index_type ghost_level)
	{
		m_ghost_level = ghost_level;
		m_node = make_unique<NodeStruct>(ni + ghost_level * 2, nj + ghost_level * 2, nk + ghost_level * 2);
		m_face = make_unique<FaceStruct>();
		m_cell = make_unique<CellStruct>();
		m_cell->Allocate(ni + ghost_level * 2 - 1, nj + ghost_level * 2 - 1, nk + ghost_level * 2 - 1);
		m_bound_map = make_unique<BoundStructManager>();
		m_idx_proxy = make_shared<IdProxyStruct>(ni + ghost_level * 2, nj + ghost_level * 2, nk + ghost_level * 2);
	}
	count_type GridStruct::GetNi()const
	{
		return m_node->GetINum();
	}

	count_type GridStruct::GetNj()const
	{
		return m_node->GetJNum();
	}

	count_type GridStruct::GetNk()const
	{
		return m_node->GetKNum();
	}
	void GridStruct::GetNodeNum(count_type& ni, count_type& nj, count_type& nk) const
	{
		ni = GetNi();
		nj = GetNj();
		nk = GetNk();
	}
	int GridStruct::GetTotalNodeNum()
	{
		return m_node->GetCount();
	}
	void GridStruct::GetRange(index_type& is, index_type& ie, index_type& js, index_type& je, index_type& ks, index_type& ke) const
	{
		is = GetGhostLevel();
		ie = GetNi() - GetGhostLevel() - 1;
		js = GetGhostLevel();
		je = GetNj() - GetGhostLevel() - 1;
		ks = GetGhostLevel();
		ke = GetNk() - GetGhostLevel() - 1;
	}
}