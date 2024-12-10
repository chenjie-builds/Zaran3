#include "GridStruct.h"
namespace zaran
{
	GridStruct::GridStruct(const string &name, Id index, Id dim) : GridBase(name, index, dim, GridType::Structured)
	{
		m_ghost_level = 1;
	}
	GridStruct::~GridStruct()
	{
	}
	void GridStruct::Allocate(Id ni, Id nj, Id nk, Id ghost_level)
	{
		m_ghost_level = ghost_level;
		m_node = std::make_shared<NodeStruct>();
		m_node->Allocate(ni + ghost_level * 2, nj + ghost_level * 2, nk + ghost_level * 2);
		m_face = std::make_shared<FaceStruct>();
		m_cell = std::make_shared<CellStruct>();
		m_cell->Allocate(ni + ghost_level * 2 - 1, nj + ghost_level * 2 - 1, nk + ghost_level * 2 - 1);
		m_bound_map = std::make_shared<BoundManagerStruct>();
		m_idx_proxy = std::make_shared<IdProxyStruct>(ni + ghost_level * 2, nj + ghost_level * 2, nk + ghost_level * 2);
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
	void GridStruct::GetNodeNum(Id&ni, Id&nj, Id&nk) const
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