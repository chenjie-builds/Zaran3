#include "GridGeneratorBlock.h"
namespace zaran
{

	void GridGeneratorBlock::CreateGrid(shared_ptr<GridBlock>& grid, GridBlockInfo& grid_info)
	{
		m_grid_info = grid_info;
		m_grid = grid;
		AllocateGrid();
		CreateNode();
		CreateBound();
		grid = m_grid;
	}
	void GridGeneratorBlock::AllocateGrid()
	{
		int ni = m_grid_info.ni;
		int nj = m_grid_info.nj;
		int nk = m_grid_info.nk;
		int ghost_size = 3;
		index_type dim = 3;
		if (nk == 1)
		{
			dim = 2;
		}
		m_grid = make_shared<GridBlock>("Structured", 1, dim);
		m_grid->Allocate(ni, nj, nk, ghost_size);
		m_grid->SetBoundBox(m_grid_info.bound_box);
	}
	void GridGeneratorBlock::CreateNode()
	{
		auto grid = GetGrid();
		auto node = grid->GetNode();
		auto ni = grid->GetNi();
		auto nj = grid->GetNj();
		auto nk = grid->GetNk();
		int ghost_size = grid->GetGhostLevel();
		auto& bound_box = grid->GetBoundBox();
		double dx = (bound_box.x_max - bound_box.x_min) / (ni - 2 * ghost_size - 1);
		double dy = (bound_box.y_max - bound_box.y_min) / (nj - 2 * ghost_size - 1);
		double dz = (bound_box.z_max - bound_box.z_min) / (nk - 2 * ghost_size - 1);
		grid->SetDx(dx);
		grid->SetDy(dy);
		grid->SetDz(dz);
		for (int idx_k = 0; idx_k < nk; ++idx_k)
		{
			for (int idx_j = 0; idx_j < nj; ++idx_j)
			{
				for (int idx_i = 0; idx_i < ni; ++idx_i)
				{
					double coord[3];
					coord[0] = bound_box.x_min + (idx_i - ghost_size) * dx;
					coord[1] = bound_box.y_min + (idx_j - ghost_size) * dy;
					coord[2] = bound_box.z_min + (idx_k - ghost_size) * dz;
					if (grid->GetDim() == 2)
						coord[2] = 0;
					node->SetCoord(idx_i, idx_j, idx_k, coord);
				}
			}
		}
	}
	void GridGeneratorBlock::CreateBound()
	{
		auto grid = GetGrid();
		auto bnd_manager = grid->GetBoundMap();
		index_type is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		// i minus
		if (!bnd_manager->HasBound(m_grid_info.bound_type_i_minus))
			bnd_manager->AllocateBound(m_grid_info.bound_type_i_minus);
		auto& bnd_im = bnd_manager->GetBound(m_grid_info.bound_type_i_minus);
		int bnd_size = (je - js + 1) * (ke - ks + 1);
		int bnd_idx = bnd_im.size();
		bnd_im.resize(bnd_im.size() + bnd_size);
		int i_bnd, j_bnd, k_bnd;
		int dir[3] = { -1, 0, 0 };
		double norm[3] = { -1, 0, 0 };
		for (int idx_k = ks; idx_k <= ke; ++idx_k)
		{
			for (int idx_j = js; idx_j <= je; ++idx_j)
			{
				i_bnd = is;
				j_bnd = idx_j;
				k_bnd = idx_k;
				auto& bnd_info = bnd_im[bnd_idx++];
				bnd_info.SetIdxSrc(i_bnd, j_bnd, k_bnd);
				bnd_info.SetDirectionSrc(dir);
				bnd_info.SetNorm(norm);
			}
		}
		// i plus
		if (!bnd_manager->HasBound(m_grid_info.bound_type_i_plus))
			bnd_manager->AllocateBound(m_grid_info.bound_type_i_plus);
		auto& bnd_ip = bnd_manager->GetBound(m_grid_info.bound_type_i_plus);
		bnd_size = (je - js + 1) * (ke - ks + 1);
		bnd_idx = bnd_ip.size();
		bnd_ip.resize(bnd_ip.size() + bnd_size);
		dir[0] = 1, dir[1] = 0, dir[2] = 0;
		norm[0] = 1, norm[1] = 0, norm[2] = 0;
		for (int idx_k = ks; idx_k <= ke; ++idx_k)
		{
			for (int idx_j = js; idx_j <= je; ++idx_j)
			{
				i_bnd = ie;
				j_bnd = idx_j;
				k_bnd = idx_k;
				auto& bnd_info = bnd_ip[bnd_idx++];
				bnd_info.SetIdxSrc(i_bnd, j_bnd, k_bnd);
				bnd_info.SetDirectionSrc(dir);
				bnd_info.SetNorm(norm);
			}
		}
		// j minus
		if (!bnd_manager->HasBound(m_grid_info.bound_type_j_minus))
			bnd_manager->AllocateBound(m_grid_info.bound_type_j_minus);
		auto& bnd_jm = bnd_manager->GetBound(m_grid_info.bound_type_j_minus);
		bnd_size = (ie - is + 1) * (ke - ks + 1);
		bnd_idx = bnd_jm.size();
		bnd_jm.resize(bnd_jm.size() + bnd_size);
		dir[0] = 0, dir[1] = -1, dir[2] = 0;
		norm[0] = 0, norm[1] = -1, norm[2] = 0;
		for (int idx_k = ks; idx_k <= ke; ++idx_k)
		{
			for (int idx_i = is; idx_i <= ie; ++idx_i)
			{
				i_bnd = idx_i;
				j_bnd = js;
				k_bnd = idx_k;
				auto& bnd_info = bnd_jm[bnd_idx++];
				bnd_info.SetIdxSrc(i_bnd, j_bnd, k_bnd);
				bnd_info.SetDirectionSrc(dir);
				bnd_info.SetNorm(norm);
			}
		}
		// j plus
		if (!bnd_manager->HasBound(m_grid_info.bound_type_j_plus))
			bnd_manager->AllocateBound(m_grid_info.bound_type_j_plus);
		auto& bnd_jp = bnd_manager->GetBound(m_grid_info.bound_type_j_plus);
		bnd_size = (ie - is + 1) * (ke - ks + 1);
		bnd_idx = bnd_jp.size();
		bnd_jp.resize(bnd_jp.size() + bnd_size);
		dir[0] = 0, dir[1] = 1, dir[2] = 0;
		norm[0] = 0, norm[1] = 1, norm[2] = 0;
		for (int idx_k = ks; idx_k <= ke; ++idx_k)
		{
			for (int idx_i = is; idx_i <= ie; ++idx_i)
			{
				i_bnd = idx_i;
				j_bnd = je;
				k_bnd = idx_k;
				auto& bnd_info = bnd_jp[bnd_idx++];
				bnd_info.SetIdxSrc(i_bnd, j_bnd, k_bnd);
				bnd_info.SetDirectionSrc(dir);
				bnd_info.SetNorm(norm);
			}
		}
		if (grid->GetDim() == 3)
		{

			// k minus
			if (!bnd_manager->HasBound(m_grid_info.bound_type_k_minus))
				bnd_manager->AllocateBound(m_grid_info.bound_type_k_minus);
			auto& bnd_km = bnd_manager->GetBound(m_grid_info.bound_type_k_minus);
			bnd_size = (ie - is + 1) * (je - js + 1);
			bnd_idx = bnd_km.size();
			bnd_km.resize(bnd_km.size() + bnd_size);
			dir[0] = 0, dir[1] = 0, dir[2] = -1;
			norm[0] = 0, norm[1] = 0, norm[2] = -1;
			for (int idx_j = js; idx_j <= je; ++idx_j)
			{
				for (int idx_i = is; idx_i <= ie; ++idx_i)
				{
					i_bnd = idx_i;
					j_bnd = idx_j;
					k_bnd = ks;
					auto& bnd_info = bnd_km[bnd_idx++];
					bnd_info.SetIdxSrc(i_bnd, j_bnd, k_bnd);
					bnd_info.SetDirectionSrc(dir);
					bnd_info.SetNorm(norm);
				}
			}
			// k plus
			if (!bnd_manager->HasBound(m_grid_info.bound_type_k_plus))
				bnd_manager->AllocateBound(m_grid_info.bound_type_k_plus);
			auto& bnd_kp = bnd_manager->GetBound(m_grid_info.bound_type_k_plus);
			bnd_size = (ie - is + 1) * (je - js + 1);
			bnd_idx = bnd_kp.size();
			bnd_kp.resize(bnd_kp.size() + bnd_size);
			dir[0] = 0, dir[1] = 0, dir[2] = 1;
			norm[0] = 0, norm[1] = 0, norm[2] = 1;
			for (int idx_j = js; idx_j <= je; ++idx_j)
			{
				for (int idx_i = is; idx_i <= ie; ++idx_i)
				{
					i_bnd = idx_i;
					j_bnd = idx_j;
					k_bnd = ke;
					auto& bnd_info = bnd_kp[bnd_idx++];
					bnd_info.SetIdxSrc(i_bnd, j_bnd, k_bnd);
					bnd_info.SetDirectionSrc(dir);
					bnd_info.SetNorm(norm);
				}
			}
		}
	}
}