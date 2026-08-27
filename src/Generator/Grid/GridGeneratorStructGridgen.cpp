#include "GridGeneratorStructGridgen.h"
#include "Log.h"
#include "ZaranError.h"
#include <fstream>
#include "MathBasic.h"
#include"Rand.h"
namespace zaran
{
	GridBuilderStructGridgen::GridBuilderStructGridgen()
	{
		// m_node_file_name = GlobalData::GetString("NodeFile");
		// m_bnd_file_name = GlobalData::GetString("BoundFile");
		std::string work_dir = GlobalData::GetString("work_dir");
		m_node_file_name = "mesh.dat";
		m_bnd_file_name = "mesh.inp";
		m_node_file_name = work_dir + "/" + m_node_file_name;
		m_bnd_file_name = work_dir + "/" + m_bnd_file_name;
	}

	void GridBuilderStructGridgen::CreateGrid(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{
		m_ghost_size = 3;
		ReadNodeFile();
		ReadBoundFile();
		SetMultiBlockInfo();
		AllocateGridMemory(grid_list);
		SetNodeCoord(grid_list);
		SetBoundInfo(grid_list);
		SetGhostNodeCoord(grid_list);
		WriteGridTest(grid_list);
	}

	void GridBuilderStructGridgen::AllocateGridMemory(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{
		grid_list.resize(GetBlockNum());
		for (int idx_block = 0; idx_block < GetBlockNum(); ++idx_block)
		{
			grid_list[idx_block] = make_shared<GridStruct>("Structured", 1, m_dim);
			auto grid = std::static_pointer_cast<GridStruct>(grid_list[idx_block]);
			grid->Allocate(m_block[idx_block].ni, m_block[idx_block].nj, m_block[idx_block].nk, m_ghost_size);
		}
	}

	void GridBuilderStructGridgen::ReadNodeFile()
	{
		std::ifstream node_file(m_node_file_name);
        //std::ofstream rand_node_file("rand_node.dat");
        //double rand_factor = 0.001;
		//随机数

		if (!node_file.is_open())
		{
			Log::error("Can't open node file:{}, Please Check!", m_node_file_name);
			throw ZaranError("Can't open node file: " + m_node_file_name);
		}
		int block_num;
		node_file >> block_num;
        //rand_node_file << block_num << std::endl;
		m_block.resize(block_num);
		for (size_t iBlock = 0; iBlock < block_num; ++iBlock)
		{
			auto& ni = m_block[iBlock].ni;
			auto& nj = m_block[iBlock].nj;
			auto& nk = m_block[iBlock].nk;
			node_file >> ni >> nj >> nk;
            //rand_node_file << ni << " " << nj << " " << nk << " " << std::endl;
			auto& coord = m_block[iBlock].node_coord;
			coord.resize(ni);
			for (int idx_i = 0; idx_i < ni; ++idx_i)
			{
				coord[idx_i].resize(nj);
				for (int idx_j = 0; idx_j < nj; ++idx_j)
				{
					coord[idx_i][idx_j].resize(nk);
				}
			}
			for (int idx_k = 0; idx_k < nk; ++idx_k)
			{
				for (int idx_j = 0; idx_j < nj; ++idx_j)
				{
					for (int idx_i = 0; idx_i < ni; ++idx_i)
					{
						node_file >> coord[idx_i][idx_j][idx_k].coord[0];
      //                  if (idx_i == 0 || idx_j == 0||idx_k==0)
						//{
      //                      rand_node_file << coord[idx_i][idx_j][idx_k].coord[0] << std::endl;
						//}
      //                  else if (idx_i == ni - 1 || idx_j == nj - 1 || idx_k == nk - 1)
						//{
						//	rand_node_file << coord[idx_i][idx_j][idx_k].coord[0] << std::endl;
						//}
						//else
						//{
      //                      double rand_x = RandDouble(-rand_factor, rand_factor);
						//	rand_node_file << coord[idx_i][idx_j][idx_k].coord[0]+rand_x << std::endl;
						//}
					}
				}
			}
			for (int idx_k = 0; idx_k < nk; ++idx_k)
			{
				for (int idx_j = 0; idx_j < nj; ++idx_j)
				{
					for (int idx_i = 0; idx_i < ni; ++idx_i)
					{
						node_file >> coord[idx_i][idx_j][idx_k].coord[1];
                        //if (idx_i == 0 || idx_j == 0 || idx_k == 0)
                        //{
                        //    rand_node_file << coord[idx_i][idx_j][idx_k].coord[1] << std::endl;
                        //}
                        //else if (idx_i == ni - 1 || idx_j == nj - 1 || idx_k == nk - 1)
                        //{
                        //    rand_node_file << coord[idx_i][idx_j][idx_k].coord[1] << std::endl;
                        //}
                        //else
                        //{
                        //    double rand_y = RandDouble(-rand_factor, rand_factor);
                        //    rand_node_file << coord[idx_i][idx_j][idx_k].coord[1] + rand_y << std::endl;
                        //}
					}
				}
			}
			for (int idx_k = 0; idx_k < nk; ++idx_k)
			{
				for (int idx_j = 0; idx_j < nj; ++idx_j)
				{
					for (int idx_i = 0; idx_i < ni; ++idx_i)
					{
						node_file >> coord[idx_i][idx_j][idx_k].coord[2];
						//if (idx_i == 0 || idx_j == 0 || idx_k == 0)
						//{
						//	rand_node_file << coord[idx_i][idx_j][idx_k].coord[2] << std::endl;
						//}
						//else if (idx_i == ni - 1 || idx_j == nj - 1 || idx_k == nk - 1)
						//{
						//	rand_node_file << coord[idx_i][idx_j][idx_k].coord[2] << std::endl;
						//}
      //                  else
      //                  {
      //                      double rand_z = RandDouble(-rand_factor, rand_factor);
      //                      rand_node_file << coord[idx_i][idx_j][idx_k].coord[2] + rand_z << std::endl;
      //                  }
					}
				}
			}
		}
		// 设置网格的维度
		if (m_block[0].nk == 1)
		{
			m_dim = 2;
		}
		else
		{
			m_dim = 3;
		}
		node_file.close();
		Log::info("Read Node File Success!");
        //rand_node_file.close();
        //exit(0);
	}
	void GridBuilderStructGridgen::ReadBoundFile()
	{
		// read inp file
		std::ifstream bnd_file(m_bnd_file_name);
		if (!bnd_file.is_open())
		{
			Log::error("Can't open bound file:{}, Please Check!", m_bnd_file_name);
			throw ZaranError("Can't open bound file: " + m_bnd_file_name);
		}
		int block_idx;
		bnd_file >> block_idx;
		int block_num;
		bnd_file >> block_num;
		for (int idx_block = 0; idx_block < block_num; ++idx_block)
		{
			int ni, nj, nk;
			bnd_file >> ni >> nj;
			if (m_dim == 3)
			{
				bnd_file >> nk;
			}
			std::string block_name;
			bnd_file >> block_name;
			int bound_num;
			bnd_file >> bound_num;
			m_block[idx_block].bound_info.resize(bound_num);
			for (int idx_bound = 0; idx_bound < bound_num; ++idx_bound)
			{
				auto& bound_info = m_block[idx_block].bound_info[idx_bound];
				bound_info.idx_block_source = idx_block;
				bnd_file >> bound_info.is_s >> bound_info.ie_s;
				bnd_file >> bound_info.js_s >> bound_info.je_s;
				if (m_dim == 3)
				{
					bnd_file >> bound_info.ks_s >> bound_info.ke_s;
				}
				else
				{
					bound_info.ks_s = 1;
					bound_info.ke_s = 1;
				}
				bound_info.dir_s[0] = bound_info.dir_s[1] = bound_info.dir_s[2] = 0;
				if (bound_info.is_s == bound_info.ie_s)
				{
					if (bound_info.is_s == 1)
						bound_info.dir_s[0] = -1;
					else
						bound_info.dir_s[0] = 1;
				}
				if (bound_info.js_s == bound_info.je_s)
				{
					if (bound_info.js_s == 1)
						bound_info.dir_s[1] = -1;
					else
						bound_info.dir_s[1] = 1;
				}
				if (m_dim == 3)
				{
					if (bound_info.ks_s == bound_info.ke_s)
					{
						if (bound_info.ks_s == 1)
							bound_info.dir_s[2] = -1;
						else
							bound_info.dir_s[2] = 1;
					}
				}
				bnd_file >> bound_info.bound_type;
				if (bound_info.bound_type == -1)
				{
					bnd_file >> bound_info.is_t >> bound_info.ie_t;
					bnd_file >> bound_info.js_t >> bound_info.je_t;
					if (m_dim == 3)
					{
						bnd_file >> bound_info.ks_t >> bound_info.ke_t;
					}
					else
					{
						bound_info.ks_t = 1;
						bound_info.ke_t = 1;
					}
					bound_info.dir_t[0] = bound_info.dir_t[1] = bound_info.dir_t[2] = 0;
					if (bound_info.is_t == bound_info.ie_t)
					{
						if (bound_info.is_t == 1)
							bound_info.dir_t[0] = -1;
						else
							bound_info.dir_t[0] = 1;
					}
					if (bound_info.js_t == bound_info.je_t)
					{
						if (bound_info.js_t == 1)
							bound_info.dir_t[1] = -1;
						else
							bound_info.dir_t[1] = 1;
					}
					if (m_dim == 3)
					{
						if (bound_info.ks_t == bound_info.ke_t)
						{
							if (bound_info.ks_t == 1)
								bound_info.dir_t[2] = -1;
							else
								bound_info.dir_t[2] = 1;
						}
					}
					bnd_file >> bound_info.block_indx_target;
					bound_info.block_indx_target -= 1;
				}
			}
		}
	}
	void GridBuilderStructGridgen::SetMultiBlockInfo()
	{
		for (int idx_block = 0; idx_block < GetBlockNum(); ++idx_block)
		{
			auto& block = m_block[idx_block];
			auto& bound = block.bound_info;
			for (int idx_bnd = 0; idx_bnd < bound.size(); ++idx_bnd)
			{
				auto& bound_info = bound[idx_bnd];
				if (bound_info.bound_type != -1)
				{
					bound_info.is_s -= 1, bound_info.ie_s -= 1;
					bound_info.js_s -= 1, bound_info.je_s -= 1;
					bound_info.ks_s -= 1, bound_info.ke_s -= 1;
					continue;
				}
				// it is i-direction boundary face
				if (bound_info.is_s == bound_info.ie_s)
				{
					if (bound_info.is_t == bound_info.ie_t)
					{
						bound_info.conn_info[0] = 0;
						if (bound_info.js_s * bound_info.js_t > 0)
						{
							bound_info.conn_info[1] = 1;
							bound_info.conn_info[2] = 2;
						}
						else
						{
							bound_info.conn_info[1] = 2;
							bound_info.conn_info[2] = 1;
						}
					}
					else if (bound_info.js_t == bound_info.je_t)
					{
						bound_info.conn_info[0] = 1;
						if (bound_info.is_s * bound_info.is_t > 0)
						{
							bound_info.conn_info[1] = 0;
							bound_info.conn_info[2] = 2;
						}
						else
						{
							bound_info.conn_info[1] = 2;
							bound_info.conn_info[2] = 0;
						}
					}
					else if (bound_info.ks_t == bound_info.ke_t)
					{
						bound_info.conn_info[0] = 2;
						if (bound_info.is_s * bound_info.is_t > 0)
						{
							bound_info.conn_info[1] = 0;
							bound_info.conn_info[2] = 1;
						}
						else
						{
							bound_info.conn_info[1] = 1;
							bound_info.conn_info[2] = 0;
						}
					}
				}
				else if (bound_info.js_s == bound_info.je_s) // it is j-direction boundary face
				{
					if (bound_info.is_t == bound_info.ie_t)
					{
						bound_info.conn_info[1] = 0;
						if (bound_info.js_s * bound_info.js_t > 0)
						{
							bound_info.conn_info[0] = 1;
							bound_info.conn_info[2] = 2;
						}
						else
						{
							bound_info.conn_info[0] = 2;
							bound_info.conn_info[2] = 1;
						}
					}
					else if (bound_info.js_t == bound_info.je_t)
					{
						bound_info.conn_info[1] = 1;
						if (bound_info.is_s * bound_info.is_t > 0)
						{
							bound_info.conn_info[0] = 0;
							bound_info.conn_info[2] = 2;
						}
						else
						{
							bound_info.conn_info[0] = 2;
							bound_info.conn_info[2] = 0;
						}
					}
					else if (bound_info.ks_t == bound_info.ke_t)
					{
						bound_info.conn_info[1] = 2;
						if (bound_info.js_s * bound_info.js_t > 0)
						{
							bound_info.conn_info[0] = 1;
							bound_info.conn_info[2] = 0;
						}
						else
						{
							bound_info.conn_info[0] = 0;
							bound_info.conn_info[2] = 1;
						}
					}
				}
				else if (bound_info.ks_s == bound_info.ke_s)
				{
					if (bound_info.is_t == bound_info.ie_t)
					{
						bound_info.conn_info[2] = 0;
						if (bound_info.is_s * bound_info.js_t > 0)
						{
							bound_info.conn_info[0] = 1;
							bound_info.conn_info[1] = 2;
						}
						else
						{
							bound_info.conn_info[0] = 2;
							bound_info.conn_info[1] = 1;
						}
					}
					else if (bound_info.js_t == bound_info.je_t)
					{
						bound_info.conn_info[2] = 1;
						if (bound_info.is_s * bound_info.is_t > 0)
						{
							bound_info.conn_info[0] = 0;
							bound_info.conn_info[1] = 2;
						}
						else
						{
							bound_info.conn_info[0] = 2;
							bound_info.conn_info[1] = 0;
						}
					}
					else if (bound_info.ks_t == bound_info.ke_t)
					{
						bound_info.conn_info[2] = 2;
						if (bound_info.is_s * bound_info.is_t > 0)
						{
							bound_info.conn_info[0] = 0;
							bound_info.conn_info[1] = 1;
						}
						else
						{
							bound_info.conn_info[0] = 1;
							bound_info.conn_info[1] = 0;
						}
					}
				}
				bound_info.is_s = abs(bound_info.is_s), bound_info.ie_s = abs(bound_info.ie_s);
				bound_info.js_s = abs(bound_info.js_s), bound_info.je_s = abs(bound_info.je_s);
				bound_info.ks_s = abs(bound_info.ks_s), bound_info.ke_s = abs(bound_info.ke_s);
				bound_info.is_t = abs(bound_info.is_t), bound_info.ie_t = abs(bound_info.ie_t);
				bound_info.js_t = abs(bound_info.js_t), bound_info.je_t = abs(bound_info.je_t);
				bound_info.ks_t = abs(bound_info.ks_t), bound_info.ke_t = abs(bound_info.ke_t);
				bound_info.is_s -= 1, bound_info.ie_s -= 1;
				bound_info.js_s -= 1, bound_info.je_s -= 1;
				bound_info.ks_s -= 1, bound_info.ke_s -= 1;
				bound_info.is_t -= 1, bound_info.ie_t -= 1;
				bound_info.js_t -= 1, bound_info.je_t -= 1;
				bound_info.ks_t -= 1, bound_info.ke_t -= 1;
				if (bound_info.is_s > bound_info.ie_s)
				{
					std::swap(bound_info.is_s, bound_info.ie_s);
					if (bound_info.conn_info[0] == 0)
					{
						std::swap(bound_info.is_t, bound_info.ie_t);
					}
					else if (bound_info.conn_info[0] == 1)
					{
						std::swap(bound_info.js_t, bound_info.je_t);
					}
					else if (bound_info.conn_info[0] == 2)
					{
						std::swap(bound_info.ks_t, bound_info.ke_t);
					}
				}
				if (bound_info.js_s > bound_info.je_s)
				{
					std::swap(bound_info.js_s, bound_info.je_s);
					if (bound_info.conn_info[0] == 0)
					{
						std::swap(bound_info.is_t, bound_info.ie_t);
					}
					else if (bound_info.conn_info[0] == 1)
					{
						std::swap(bound_info.js_t, bound_info.je_t);
					}
					else if (bound_info.conn_info[0] == 2)
					{
						std::swap(bound_info.ks_t, bound_info.ke_t);
					}
				}
				if (bound_info.ks_s > bound_info.ke_s)
				{
					std::swap(bound_info.ks_s, bound_info.ke_s);
					if (bound_info.conn_info[0] == 0)
					{
						std::swap(bound_info.is_t, bound_info.ie_t);
					}
					else if (bound_info.conn_info[0] == 1)
					{
						std::swap(bound_info.js_t, bound_info.je_t);
					}
					else if (bound_info.conn_info[0] == 2)
					{
						std::swap(bound_info.ks_t, bound_info.ke_t);
					}
				}
			}
		}
	}
	void GridBuilderStructGridgen::SetBoundInfo(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{
		std::map<int, string> gridgen_bound = {
			{-1, "connection"},
			{0, "none"},
			{1, "interblock_connection"},
			{2, "wall"},
			{3, "symmetry"},
			{4, "farfield"},
			{5, "inlet"},
			{6, "outlet"} };
		for (int block_indx = 0; block_indx < GetBlockNum(); ++block_indx)
		{
			auto grid_src = std::static_pointer_cast<GridStruct>(grid_list[block_indx]);
			auto bnd_manager = grid_src->GetBoundMap();
			auto node_src = grid_src->GetNode();
			int i_bnd_src, j_bnd_src, k_bns_src;
			int i_bnd_tgt, j_bnd_tgt, k_bns_tgt;
			double norm[3];
			int iter_bnd_s = 0, iter_bnd_e = 0;
			for (auto& bound_info : m_block[block_indx].bound_info)
			{
				for (int k = bound_info.ks_s; k <= bound_info.ke_s; ++k)
				{
					for (int j = bound_info.js_s; j <= bound_info.je_s; ++j)
					{
						for (int i = bound_info.is_s; i <= bound_info.ie_s; ++i)
						{
							i_bnd_src = i + m_ghost_size;
							j_bnd_src = j + m_ghost_size;
							k_bns_src = k + m_ghost_size;
							auto ref_node_coord = node_src->GetCoord(i_bnd_src - bound_info.dir_s[0], j_bnd_src - bound_info.dir_s[1], k_bns_src - bound_info.dir_s[2]);
							auto bound_node_coord = node_src->GetCoord(i_bnd_src, j_bnd_src, k_bns_src);
							// 根据同层的邻居节点计算法向量
							dynamic_array<const double*> neighbor_node_coord;
							neighbor_node_coord.push_back(node_src->GetCoord(i_bnd_src, j_bnd_src, k_bns_src));
							if (bound_info.ks_s != bound_info.ke_s)
							{
								if (k != bound_info.ks_s)
								{
									neighbor_node_coord.push_back(node_src->GetCoord(i_bnd_src, j_bnd_src, k_bns_src - 1));
								}
								if (k != bound_info.ke_s)
								{
									neighbor_node_coord.push_back(node_src->GetCoord(i_bnd_src, j_bnd_src, k_bns_src + 1));
								}
							}
							if (bound_info.js_s != bound_info.je_s)
							{
								if (j != bound_info.js_s)
								{
									neighbor_node_coord.push_back(node_src->GetCoord(i_bnd_src, j_bnd_src - 1, k_bns_src));
								}
								if (j != bound_info.je_s)
								{
									neighbor_node_coord.push_back(node_src->GetCoord(i_bnd_src, j_bnd_src + 1, k_bns_src));
								}
							}
							if (bound_info.is_s != bound_info.ie_s)
							{
								if (i != bound_info.is_s)
								{
									neighbor_node_coord.push_back(node_src->GetCoord(i_bnd_src - 1, j_bnd_src, k_bns_src));
								}
								if (i != bound_info.ie_s)
								{
									neighbor_node_coord.push_back(node_src->GetCoord(i_bnd_src + 1, j_bnd_src, k_bns_src));
								}
							}
							if (grid_src->GetDim() == TWO_DIM)
							{
								double a, b, c;
								LineFit2D(neighbor_node_coord.data(), neighbor_node_coord.size(), a, b, c);
								norm[0] = a;
								norm[1] = b;
								norm[2] = 0;
							}
							else
							{
								double a, b, c, d;
								PlaneFit3D(neighbor_node_coord.data(), neighbor_node_coord.size(), a, b, c, d);
								norm[0] = a;
								norm[1] = b;
								norm[2] = c;
							}
							// 边界法向量方向与参考节点到边界节点方向相反,保证法向量指向外部
							double bnd_to_ref[3];
							bnd_to_ref[0] = ref_node_coord[0] - bound_node_coord[0];
							bnd_to_ref[1] = ref_node_coord[1] - bound_node_coord[1];
							bnd_to_ref[2] = ref_node_coord[2] - bound_node_coord[2];
							double norm_dot_bnd_to_ref = norm[0] * bnd_to_ref[0] + norm[1] * bnd_to_ref[1] + norm[2] * bnd_to_ref[2];
							if (norm_dot_bnd_to_ref > 0)
							{
								norm[0] = -norm[0];
								norm[1] = -norm[1];
								norm[2] = -norm[2];
							}
							double norm_len = sqrt(norm[0] * norm[0] + norm[1] * norm[1] + norm[2] * norm[2]);
							norm[0] /= norm_len;
							norm[1] /= norm_len;
							norm[2] /= norm_len;

							if (bound_info.bound_type < 0)
							{
								int iter_i, iter_j, iter_k;
								if (bound_info.conn_info[0] == 0)
								{
									iter_i = i - bound_info.is_s;
								}
								else if (bound_info.conn_info[0] == 1)
								{
									iter_j = i - bound_info.is_s;
								}
								else if (bound_info.conn_info[0] == 2)
								{
									iter_k = i - bound_info.is_s;
								}
								if (bound_info.conn_info[1] == 0)
								{
									iter_i = j - bound_info.js_s;
								}
								else if (bound_info.conn_info[1] == 1)
								{
									iter_j = j - bound_info.js_s;
								}
								else if (bound_info.conn_info[1] == 2)
								{
									iter_k = j - bound_info.js_s;
								}
								if (bound_info.conn_info[2] == 0)
								{
									iter_i = k - bound_info.ks_s;
								}
								else if (bound_info.conn_info[2] == 1)
								{
									iter_j = k - bound_info.ks_s;
								}
								else if (bound_info.conn_info[2] == 2)
								{
									iter_k = k - bound_info.ks_s;
								}
								i_bnd_tgt = bound_info.is_t + iter_i * Sign(bound_info.ie_t - bound_info.is_t) + m_ghost_size;
								j_bnd_tgt = bound_info.js_t + iter_j * Sign(bound_info.je_t - bound_info.js_t) + m_ghost_size;
								k_bns_tgt = bound_info.ks_t + iter_k * Sign(bound_info.ke_t - bound_info.ks_t) + m_ghost_size;
								auto grid_tgt = std::static_pointer_cast<GridStruct>(grid_list[bound_info.block_indx_target]);
								auto node_tgt = grid_tgt->GetNode();
								for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
								{
									i_bnd_src += bound_info.dir_s[0];
									j_bnd_src += bound_info.dir_s[1];
									k_bns_src += bound_info.dir_s[2];
									i_bnd_tgt -= bound_info.dir_t[0];
									j_bnd_tgt -= bound_info.dir_t[1];
									k_bns_tgt -= bound_info.dir_t[2];
									BoundStruct bound(i_bnd_src, j_bnd_src, k_bns_src, bound_info.dir_s, norm, i_bnd_tgt, j_bnd_tgt, k_bns_tgt, bound_info.dir_t, bound_info.block_indx_target);
									bnd_manager->AddBound(gridgen_bound[bound_info.bound_type], bound);
									node_src->SetCoord(i_bnd_src, j_bnd_src, k_bns_src, node_tgt->GetCoord(i_bnd_tgt, j_bnd_tgt, k_bns_tgt));

								}
								BoundStruct bound(i_bnd_src, j_bnd_src, k_bns_src, bound_info.dir_s, norm, i_bnd_tgt, j_bnd_tgt, k_bns_tgt, bound_info.dir_t, bound_info.block_indx_target);
								bnd_manager->AddBound(gridgen_bound[bound_info.bound_type], bound);

							}
							else
							{
								BoundStruct bound(i_bnd_src, j_bnd_src, k_bns_src, bound_info.dir_s, norm);
								bnd_manager->AddBound(gridgen_bound[bound_info.bound_type], bound);
							}
						}
					}
				}
			}
		}
	}
	void GridBuilderStructGridgen::SetNodeCoord(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{
		for (int idx_block = 0; idx_block < GetBlockNum(); ++idx_block)
		{
			auto grid = std::static_pointer_cast<GridStruct>(grid_list[idx_block]);
			auto node = grid->GetNode();
			index_type is, ie, js, je, ks, ke;
			grid->GetRange(is, ie, js, je, ks, ke);
			for (int idx_k = ks; idx_k <= ke; ++idx_k)
			{
				for (int idx_j = js; idx_j <= je; ++idx_j)
				{
					for (int idx_i = is; idx_i <= ie; ++idx_i)
					{
						node->SetCoord(idx_i, idx_j, idx_k, m_block[idx_block].node_coord[idx_i - m_ghost_size][idx_j - m_ghost_size][idx_k - m_ghost_size].coord);
					}
				}
			}
		}
	}
	void GridBuilderStructGridgen::SetGhostNodeCoord3D(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{
		for (int block_indx = 0; block_indx < GetBlockNum(); ++block_indx)
		{
			auto grid_src = std::static_pointer_cast<GridStruct>(grid_list[block_indx]);
			auto bnd_manager = grid_src->GetBoundMap();
			auto node_src = grid_src->GetNode();
			int i_bnd_src, j_bnd_src, k_bns_src;
			int i_bnd_tgt, j_bnd_tgt, k_bns_tgt;
			int i_ghost, j_ghost, k_ghost;
			double norm[3];
			int iter_bnd_s, iter_bnd_e;
			const double* bnd_node_coord, * ref_node_coord;
			double ghost_node_coord[3];
			index_type is, ie, js, je, ks, ke;
			grid_src->GetRange(is, ie, js, je, ks, ke);
			// step 1: set boundary ghost node coord
			for (auto& bound_info : m_block[block_indx].bound_info)
			{
				Log::info("Bound Type: {}", bound_info.bound_type);
				for (int k = bound_info.ks_s; k <= bound_info.ke_s; ++k)
				{
					for (int j = bound_info.js_s; j <= bound_info.je_s; ++j)
					{
						for (int i = bound_info.is_s; i <= bound_info.ie_s; ++i)
						{
							i_bnd_src = i + m_ghost_size;
							j_bnd_src = j + m_ghost_size;
							k_bns_src = k + m_ghost_size;

							auto grid_tgt = std::static_pointer_cast<GridStruct>(grid_list[bound_info.block_indx_target]);
							auto node_tgt = grid_tgt->GetNode();
							for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
							{
								i_ghost = i_bnd_src + (iGhost + 1) * bound_info.dir_s[0];
								j_ghost = j_bnd_src + (iGhost + 1) * bound_info.dir_s[1];
								k_ghost = k_bns_src + (iGhost + 1) * bound_info.dir_s[2];
								if (bound_info.bound_type == -1)
								{
									int iter_i, iter_j, iter_k;
									if (bound_info.conn_info[0] == 0)
									{
										iter_i = i - bound_info.is_s;
									}
									else if (bound_info.conn_info[0] == 1)
									{
										iter_j = i - bound_info.is_s;
									}
									else if (bound_info.conn_info[0] == 2)
									{
										iter_k = i - bound_info.is_s;
									}
									if (bound_info.conn_info[1] == 0)
									{
										iter_i = j - bound_info.js_s;
									}
									else if (bound_info.conn_info[1] == 1)
									{
										iter_j = j - bound_info.js_s;
									}
									else if (bound_info.conn_info[1] == 2)
									{
										iter_k = j - bound_info.js_s;
									}
									if (bound_info.conn_info[2] == 0)
									{
										iter_i = k - bound_info.ks_s;
									}
									else if (bound_info.conn_info[2] == 1)
									{
										iter_j = k - bound_info.ks_s;
									}
									else if (bound_info.conn_info[2] == 2)
									{
										iter_k = k - bound_info.ks_s;
									}
									i_bnd_tgt = bound_info.is_t + iter_i * Sign(bound_info.ie_t - bound_info.is_t) + m_ghost_size;
									j_bnd_tgt = bound_info.js_t + iter_j * Sign(bound_info.je_t - bound_info.js_t) + m_ghost_size;
									k_bns_tgt = bound_info.ks_t + iter_k * Sign(bound_info.ke_t - bound_info.ks_t) + m_ghost_size;
									auto src_bnd_coord = node_src->GetCoord(i_bnd_src, j_bnd_src, k_bns_src);
									auto tgt_bnd_coord = node_tgt->GetCoord(i_bnd_tgt, j_bnd_tgt, k_bns_tgt);
									ref_node_coord = node_tgt->GetCoord(i_bnd_tgt - bound_info.dir_t[0] * (iGhost + 1), j_bnd_tgt - bound_info.dir_t[1] * (iGhost + 1), k_bns_tgt - bound_info.dir_t[2] * (iGhost + 1));
									for (int i = 0; i < 3; ++i)
									{
										ghost_node_coord[i] = ref_node_coord[i];
									}
								}
								else
								{
									ref_node_coord = node_src->GetCoord(i_bnd_src - bound_info.dir_s[0] * (iGhost + 1), j_bnd_src - bound_info.dir_s[1] * (iGhost + 1), k_bns_src - bound_info.dir_s[2] * (iGhost + 1));
									for (int i = 0; i < 3; ++i)
									{
										ghost_node_coord[i] = 2 * node_src->GetCoord(i_bnd_src, j_bnd_src, k_bns_src)[i] - ref_node_coord[i];
									}
								}
								node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
							}
						}
					}
				}
			}
			// step 2: set edge ghost node coord
			for (int i = is; i <= ie; ++i)
			{
				i_ghost = i;
				for (int jGhost = 0; jGhost < m_ghost_size; ++jGhost)
				{
					for (int kGhost = 0; kGhost < m_ghost_size; ++kGhost)
					{
						j_ghost = m_ghost_size - (jGhost + 1);
						k_ghost = m_ghost_size - (kGhost + 1);
						double ghost_node_coord[3];
						ghost_node_coord[0] = 0.5 * (2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[0] -
							node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[0] +
							2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[0] -
							node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[0]);
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] -
							node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[2] -
							node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

						j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
						k_ghost = m_ghost_size - (kGhost + 1);
						ghost_node_coord[0] = 0.5 * (2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[0] -
							node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[0] +
							2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[0] -
							node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[0]);
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] -
							node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[2] -
							node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

						j_ghost = m_ghost_size - (jGhost + 1);
						k_ghost = (grid_src->GetNk() - m_ghost_size - 1) + (kGhost + 1);
						ghost_node_coord[0] = 0.5 * (2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[0] -
							node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[0] +
							2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[0] -
							node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[0]);
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] -
							node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[2] -
							node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

						j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
						k_ghost = (grid_src->GetNk() - m_ghost_size - 1) + (kGhost + 1);
						ghost_node_coord[0] = 0.5 * (2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[0] -
							node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[0] +
							2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[0] -
							node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[0]);
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] -
							node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[2] -
							node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
					}
				}
			}
			for (int j = js; j <= je; ++j)
			{
				j_ghost = j;
				for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
				{
					for (int kGhost = 0; kGhost < m_ghost_size; ++kGhost)
					{
						i_ghost = m_ghost_size - (iGhost + 1);
						k_ghost = m_ghost_size - (kGhost + 1);
						double ghost_node_coord[3];
						ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] -
							node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 0.5 * (2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[1] -
							node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[1] +
							2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[1] -
							node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[1]);
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[2] -
							node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
						k_ghost = m_ghost_size - (kGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] -
							node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 0.5 * (2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[1] -
							node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[1] +
							2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[1] -
							node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[1]);
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[2] -
							node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						i_ghost = m_ghost_size - (iGhost + 1);
						k_ghost = (grid_src->GetNk() - m_ghost_size - 1) + (kGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] -
							node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 0.5 * (2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[1] -
							node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[1] +
							2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[1] -
							node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[1]);
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[2] -
							node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
						k_ghost = (grid_src->GetNk() - m_ghost_size - 1) + (kGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] -
							node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 0.5 * (2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[1] -
							node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[1] +
							2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[1] -
							node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[1]);
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[2] -
							node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
					}
				}
			}
			for (int k = ks; k <= ke; ++k)
			{
				k_ghost = k;
				for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
				{
					for (int jGhost = 0; jGhost < m_ghost_size; ++jGhost)
					{
						i_ghost = m_ghost_size - (iGhost + 1);
						j_ghost = m_ghost_size - (jGhost + 1);
						if (i_ghost == 2 && j_ghost == 0 && k_ghost == 4)
							double ghost_node_coord[3];
						ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] -
							node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] -
							node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 0.5 * (2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[2] -
							node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[2] +
							2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[2] -
							node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[2]);
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

						i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
						j_ghost = m_ghost_size - (jGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] -
							node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] -
							node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 0.5 * (2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[2] -
							node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[2] +
							2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[2] -
							node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[2]);
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

						i_ghost = m_ghost_size - (iGhost + 1);
						j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] -
							node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] -
							node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 0.5 * (2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[2] -
							node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[2] +
							2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[2] -
							node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[2]);
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

						i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
						j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] -
							node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] -
							node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 0.5 * (2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[2] -
							node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[2] +
							2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[2] -
							node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[2]);
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
					}
				}
			}

			// step 3: set corner ghost node coord
			for (int kGhost = 0; kGhost < m_ghost_size; ++kGhost)
			{
				for (int jGhost = 0; jGhost < m_ghost_size; ++jGhost)
				{
					for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
					{
						i_ghost = m_ghost_size - (iGhost + 1);
						j_ghost = m_ghost_size - (jGhost + 1);
						k_ghost = m_ghost_size - (kGhost + 1);
						double ghost_node_coord[3];
						ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] - node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] - node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[2] - node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

						i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
						j_ghost = m_ghost_size - (jGhost + 1);
						k_ghost = m_ghost_size - (kGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] - node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] - node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[2] - node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						i_ghost = m_ghost_size - (iGhost + 1);
						j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
						k_ghost = m_ghost_size - (kGhost + 1);
						if (i_ghost == 2 && j_ghost == 2 && k_ghost == 16)
							Log::info("----");
						ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] - node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] - node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[2] - node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
						j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
						k_ghost = m_ghost_size - (kGhost + 1);
						if (i_ghost == 2 && j_ghost == 2 && k_ghost == 16)
							Log::info("----");
						ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] - node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] - node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, m_ghost_size)[2] - node_src->GetCoord(i_ghost, j_ghost, m_ghost_size + (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						i_ghost = m_ghost_size - (iGhost + 1);
						j_ghost = m_ghost_size - (jGhost + 1);
						k_ghost = (grid_src->GetNk() - m_ghost_size - 1) + (kGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] - node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] - node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[2] - node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
						j_ghost = m_ghost_size - (jGhost + 1);
						k_ghost = (grid_src->GetNk() - m_ghost_size - 1) + (kGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] - node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] - node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[2] - node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						i_ghost = m_ghost_size - (iGhost + 1);
						j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
						k_ghost = (grid_src->GetNk() - m_ghost_size - 1) + (kGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] - node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] - node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[2] - node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
						j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
						k_ghost = (grid_src->GetNk() - m_ghost_size - 1) + (kGhost + 1);
						ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] - node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
						ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] - node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
						ghost_node_coord[2] = 2 * node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1)[2] - node_src->GetCoord(i_ghost, j_ghost, grid_src->GetNk() - m_ghost_size - 1 - (kGhost + 1))[2];
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
					}
				}
			}
		}
	}

	void GridBuilderStructGridgen::SetGhostNodeCoord(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{
		if (m_dim == 2)
		{
			SetGhostNodeCoord2D(grid_list);
		}
		else if (m_dim == 3)
		{
			SetGhostNodeCoord3D(grid_list);
		}
	}

	void GridBuilderStructGridgen::SetGhostNodeCoord2D(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{
		// step 1: set boundary ghost node coord
		for (int block_indx = 0; block_indx < GetBlockNum(); ++block_indx)
		{
			auto grid_src = std::static_pointer_cast<GridStruct>(grid_list[block_indx]);

			auto bound_map = grid_src->GetBoundMap();
			auto node_src = grid_src->GetNode();
			index_type i_bnd_src, j_bnd_src, k_bnd_src;
			index_type i_bnd_tgt, j_bnd_tgt, k_bns_tgt;
			index_type i_ref, j_ref, k_ref;
			index_type i_ghost, j_ghost, k_ghost;
			double norm[3];
			int iter_bnd_s, iter_bnd_e;
			const double* bnd_node_coord, * ref_node_coord;
			double ghost_node_coord[3];
			for (auto& boundary : bound_map->GetBoundMap())
			{
				auto& bound_name = boundary.first;
				auto& bound = boundary.second;
				Log::info("Bound Name: {}", bound_name);
				if (bound_name == "connection")
				{
					for (int iBound = 0; iBound < bound.size(); ++iBound)
					{
						bound[iBound].GetIdx(i_bnd_src, j_bnd_src, k_bnd_src);
						for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
						{
							auto dir_src = bound[iBound].GetDirectionSrc();
							i_ghost = i_bnd_src + dir_src[0] * (iGhost + 1);
							j_ghost = i_bnd_src + dir_src[1] * (iGhost + 1);
							k_ghost = i_bnd_src + dir_src[2] * (iGhost + 1);
							auto dir_tgt = bound[iBound].GetDirectionTgt();
							auto src_bnd_coord = node_src->GetCoord(i_bnd_src, j_bnd_src, k_bnd_src);
							bound[iBound].GetIdxTgt(i_bnd_tgt, j_bnd_tgt, k_bns_tgt);
							auto grid_tgt = std::static_pointer_cast<GridStruct>(grid_list[bound[iBound].GetTargetBlock()]);
							auto node_tgt = grid_tgt->GetNode();
							auto tgt_bnd_coord = node_tgt->GetCoord(i_bnd_tgt, j_bnd_tgt, k_bns_tgt);
							i_ref = i_bnd_tgt - dir_tgt[0] * (iGhost + 1);
							j_ref = j_bnd_tgt - dir_tgt[1] * (iGhost + 1);
							k_ref = k_bns_tgt - dir_tgt[2] * (iGhost + 1);
							ref_node_coord = node_tgt->GetCoord(i_ref, j_ref, k_ref);
							for (int i = 0; i < 3; ++i)
							{
								ghost_node_coord[i] = ref_node_coord[i];
							}
							node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

						}
					}
				}
				else
				{
					for (int iBound = 0; iBound < bound.size(); ++iBound)
					{
						bound[iBound].GetIdx(i_bnd_src, j_bnd_src, k_bnd_src);
						auto bound_norm = bound[iBound].GetNorm();
						auto dir_src = bound[iBound].GetDirectionSrc();
						for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
						{
							i_ghost = i_bnd_src + dir_src[0] * (iGhost + 1);
							j_ghost = j_bnd_src + dir_src[1] * (iGhost + 1);
							k_ghost = k_bnd_src + dir_src[2] * (iGhost + 1);
							i_ref = i_bnd_src - dir_src[0] * (iGhost + 1);
							j_ref = j_bnd_src - dir_src[1] * (iGhost + 1);
							k_ref = k_bnd_src - dir_src[2] * (iGhost + 1);
							bnd_node_coord = node_src->GetCoord(i_bnd_src, j_bnd_src, k_bnd_src);
							ref_node_coord = node_src->GetCoord(i_ref, j_ref, k_ref);
							double vec_bnd_ref[3];
							double proj_norm = 0.0;
							for (int i = 0; i < 3; ++i)
							{
								vec_bnd_ref[i] = ref_node_coord[i] - bnd_node_coord[i];
								proj_norm += vec_bnd_ref[i] * bound_norm[i];
							}
							for (int i = 0; i < 3; ++i)
							{
								//ghost_node_coord[i] = bnd_node_coord[i] + (vec_bnd_ref[i] - 2.0 * proj_norm * bound_norm[i]);
								ghost_node_coord[i] = bnd_node_coord[i] - (vec_bnd_ref[i] );
							}
							node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						}
					}
				}
			}
		}
		//return;



		for (int block_indx = 0; block_indx < GetBlockNum(); ++block_indx)
		{
			auto grid_src = std::static_pointer_cast<GridStruct>(grid_list[block_indx]);
			auto bnd_manager = grid_src->GetBoundMap();
			auto node_src = grid_src->GetNode();
			int i_bnd_src, j_bnd_src, k_bnd_src;
			int i_bnd_tgt, j_bnd_tgt, k_bnd_tgt;
			int i_ghost, j_ghost, k_ghost;
			double norm[3];
			int iter_bnd_s, iter_bnd_e;
			const double* bnd_node_coord, * ref_node_coord;
			double ghost_node_coord[3];
			 //step 1: set boundary ghost node coord
			/*for (auto& bound_info : m_block[block_indx].bound_info)
			{
				Log::info("Bound Type: {}", bound_info.bound_type);
				for (int k = bound_info.ks_s; k <= bound_info.ke_s; ++k)
				{
					for (int j = bound_info.js_s; j <= bound_info.je_s; ++j)
					{
						for (int i = bound_info.is_s; i <= bound_info.ie_s; ++i)
						{
							i_bnd_src = i + m_ghost_size;
							j_bnd_src = j + m_ghost_size;
							k_bnd_src = k + m_ghost_size;

							auto grid_tgt = dynamic_cast<GridStruct*>(grid_list[bound_info.block_indx_target]);
							auto node_tgt = grid_tgt->GetNode();
							for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
							{
								i_ghost = i_bnd_src + (iGhost + 1) * bound_info.dir_s[0];
								j_ghost = j_bnd_src + (iGhost + 1) * bound_info.dir_s[1];
								k_ghost = k_bnd_src + (iGhost + 1) * bound_info.dir_s[2];
								if (bound_info.bound_type == -1)
								{
									int iter_i, iter_j, iter_k;
									if (bound_info.conn_info[0] == 0)
									{
										iter_i = i - bound_info.is_s;
									}
									else if (bound_info.conn_info[0] == 1)
									{
										iter_j = i - bound_info.is_s;
									}
									else if (bound_info.conn_info[0] == 2)
									{
										iter_k = i - bound_info.is_s;
									}
									if (bound_info.conn_info[1] == 0)
									{
										iter_i = j - bound_info.js_s;
									}
									else if (bound_info.conn_info[1] == 1)
									{
										iter_j = j - bound_info.js_s;
									}
									else if (bound_info.conn_info[1] == 2)
									{
										iter_k = j - bound_info.js_s;
									}
									if (bound_info.conn_info[2] == 0)
									{
										iter_i = k - bound_info.ks_s;
									}
									else if (bound_info.conn_info[2] == 1)
									{
										iter_j = k - bound_info.ks_s;
									}
									else if (bound_info.conn_info[2] == 2)
									{
										iter_k = k - bound_info.ks_s;
									}
									i_bnd_tgt = bound_info.is_t + iter_i * Sign(bound_info.ie_t - bound_info.is_t) + m_ghost_size;
									j_bnd_tgt = bound_info.js_t + iter_j * Sign(bound_info.je_t - bound_info.js_t) + m_ghost_size;
									k_bnd_tgt = bound_info.ks_t + iter_k * Sign(bound_info.ke_t - bound_info.ks_t) + m_ghost_size;
									auto src_bnd_coord = node_src->GetCoord(i_bnd_src, j_bnd_src, k_bnd_src);
									auto tgt_bnd_coord = node_tgt->GetCoord(i_bnd_tgt, j_bnd_tgt, k_bnd_tgt);
									ref_node_coord = node_tgt->GetCoord(i_bnd_tgt - bound_info.dir_t[0] * (iGhost + 1), j_bnd_tgt - bound_info.dir_t[1] * (iGhost + 1), k_bnd_tgt - bound_info.dir_t[2] * (iGhost + 1));
									for (int i = 0; i < 3; ++i)
									{
										ghost_node_coord[i] = ref_node_coord[i];
									}
								}
								else
								{
									ref_node_coord = node_src->GetCoord(i_bnd_src - bound_info.dir_s[0] * (iGhost + 1), j_bnd_src - bound_info.dir_s[1] * (iGhost + 1), k_bnd_src - bound_info.dir_s[2] * (iGhost + 1));
									for (int i = 0; i < 3; ++i)
									{
										ghost_node_coord[i] = 2 * node_src->GetCoord(i_bnd_src, j_bnd_src, k_bnd_src)[i] - ref_node_coord[i];
									}
								}
								node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
							}
						}
					}
				}
			}*/
			// step 2: set corner ghost node coord
			for (int jGhost = 0; jGhost < m_ghost_size; ++jGhost)
			{
				for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
				{
					i_ghost = m_ghost_size - (iGhost + 1);
					j_ghost = m_ghost_size - (jGhost + 1);
					k_ghost = m_ghost_size;
					double ghost_node_coord[3];
					ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] - node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
					ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] - node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
					ghost_node_coord[2] = node_src->GetCoord(m_ghost_size, m_ghost_size, m_ghost_size)[2];
					node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

					i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
					j_ghost = m_ghost_size - (jGhost + 1);
					k_ghost = m_ghost_size;
					ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] - node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
					ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, m_ghost_size, k_ghost)[1] - node_src->GetCoord(i_ghost, m_ghost_size + (jGhost + 1), k_ghost)[1];
					ghost_node_coord[2] = node_src->GetCoord(m_ghost_size, m_ghost_size, m_ghost_size)[2];
					node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

					i_ghost = m_ghost_size - (iGhost + 1);
					j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
					k_ghost = m_ghost_size;
					ghost_node_coord[0] = 2 * node_src->GetCoord(m_ghost_size, j_ghost, k_ghost)[0] - node_src->GetCoord(m_ghost_size + (iGhost + 1), j_ghost, k_ghost)[0];
					ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] - node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
					ghost_node_coord[2] = node_src->GetCoord(m_ghost_size, m_ghost_size, m_ghost_size)[2];
					node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);

					i_ghost = (grid_src->GetNi() - m_ghost_size - 1) + (iGhost + 1);
					j_ghost = (grid_src->GetNj() - m_ghost_size - 1) + (jGhost + 1);
					k_ghost = m_ghost_size;
					ghost_node_coord[0] = 2 * node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1, j_ghost, k_ghost)[0] - node_src->GetCoord(grid_src->GetNi() - m_ghost_size - 1 - (iGhost + 1), j_ghost, k_ghost)[0];
					ghost_node_coord[1] = 2 * node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1, k_ghost)[1] - node_src->GetCoord(i_ghost, grid_src->GetNj() - m_ghost_size - 1 - (jGhost + 1), k_ghost)[1];
					ghost_node_coord[2] = node_src->GetCoord(m_ghost_size, m_ghost_size, m_ghost_size)[2];
					node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
				}
			}

			int ni, nj;
			ni = grid_src->GetNi();
			nj = grid_src->GetNj();
			for (int i = 0; i < ni; ++i)
			{
				for (int j = 0; j < nj; ++j)
				{
					i_bnd_src = i;
					j_bnd_src = j;
					k_bnd_src = m_ghost_size;
					for (int iGhost = 0; iGhost < m_ghost_size; ++iGhost)
					{
						i_ghost = i_bnd_src;
						j_ghost = j_bnd_src;
						k_ghost = k_bnd_src + (iGhost + 1);
						ghost_node_coord[0] = node_src->GetCoord(i_bnd_src, j_bnd_src, k_bnd_src)[0];
						ghost_node_coord[1] = node_src->GetCoord(i_bnd_src, j_bnd_src, k_bnd_src)[1];
						ghost_node_coord[2] = node_src->GetCoord(i_bnd_src, j_bnd_src, k_bnd_src)[2] + 1.0 * (iGhost + 1);
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
						k_ghost = k_bnd_src - (iGhost + 1);
						ghost_node_coord[2] = node_src->GetCoord(i_bnd_src, j_bnd_src, k_bnd_src)[2] - 1.0 * (iGhost + 1);
						node_src->SetCoord(i_ghost, j_ghost, k_ghost, ghost_node_coord);
					}
				}
			}
		}
	}

	void GridBuilderStructGridgen::WriteGridTest(dynamic_array<shared_ptr<GridBase>>& grid_list)
	{

		for (int idx_block = 0; idx_block < GetBlockNum(); ++idx_block)
		{
			std::ofstream file("grid_test" + std::to_string(idx_block) + ".dat");
			auto grid = std::static_pointer_cast<GridStruct>(grid_list[idx_block]);
			auto node = grid->GetNode();
			int ni, nj, nk;
			ni = grid->GetNi();
			nj = grid->GetNj();
			nk = grid->GetNk();
			index_type is, ie, js, je, ks, ke;
			grid->GetRange(is, ie, js, je, ks, ke);
			file << "TITLE=\"Flow Field\"\n";
			file << "VARIABLES=\"X\",\"Y\",\"Z\",\"I\",\"J\",\"K\"\n";
			file << "ZONE I=" << ni << ", J=" << nj << ", K=" << nk << ", F=POINT\n";
			for (int k = 0; k < nk; ++k)
			{
				for (int j = 0; j < nj; ++j)
				{
					for (int i = 0; i < ni; ++i)
					{
						auto coord = node->GetCoord(i, j, k);
						if (abs(coord[0]) > 100)
							file << 0 << " ";
						else
							file << coord[0] << " ";
						if (abs(coord[1]) > 100)

							file << 0 << " ";
						else
							file << coord[1] << " ";
						if (abs(coord[2]) > 100)
							file << 0 << " ";
						else
							file << coord[2] << " ";
						file << " " << i << " " << j << " " << k << std::endl;
					}
				}
			}
		}
	}

} // namespace zaran