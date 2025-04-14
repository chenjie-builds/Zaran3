#include "GridConvert.h"
namespace zaran
{


	void GridConvertPiflow2Flexible::ReadFile()
	{
		std::ifstream file(m_file_name);
		if (!file.is_open())
		{
			throw std::runtime_error("file not found");
		}

		ReadNodeCoord(file);
		ReadInnerNode(file);
		ReadGradNode(file);
		ReadBound(file);
		ReadCell(file);

	}

	void GridConvertPiflow2Flexible::WriteFile()
	{
		string node_file = "node.dat";
		string cell_file = "cell.dat";

		std::ofstream file(node_file);
		if (!file.is_open())
		{
			throw std::runtime_error("file not found");
		}
		WriteNodeCoord(file);
		WriteInnerNode(file);
		WriteBound(file);
		file.close();
		file.open(cell_file);
		if (!file.is_open())
		{
			throw std::runtime_error("file not found");
		}
		WriteCell(file);
		file.close();
	}


	void GridConvertPiflow2Flexible::ReadNodeCoord(std::ifstream& file)
	{
		index_type node_count;
		file >> node_count;
		m_coord.resize(node_count);
		for (index_type i = 0; i < node_count; i++)
		{
			file >> m_coord[i].x >> m_coord[i].y >> m_coord[i].z;
		}

	}

	void GridConvertPiflow2Flexible::ReadInnerNode(std::ifstream& file)
	{
		index_type inner_node_count;
		file >> inner_node_count;
		m_inner_node.resize(inner_node_count);
		for (index_type i = 0; i < inner_node_count; i++)
		{
			file >> m_inner_node[i].idx;
			for (index_type j = 0; j < 6; j++)
			{
				file >> m_inner_node[i].neighbor[j];
			}
		}

	}

	void GridConvertPiflow2Flexible::ReadGradNode(std::ifstream& file)
	{
		index_type grad_node_count;
		file >> grad_node_count;
		m_grad_node.resize(grad_node_count);
		for (index_type i = 0; i < grad_node_count; i++)
		{
			file >> m_grad_node[i].idx;
			index_type neighbor_count;
			file >> neighbor_count;
			m_grad_node[i].neighbor.resize(neighbor_count);
			for (index_type j = 0; j < neighbor_count; j++)
			{
				file >> m_grad_node[i].neighbor[j];
			}
		}

	}

	void GridConvertPiflow2Flexible::ReadBound(std::ifstream& file)
	{
		index_type bound_count;
		file >> bound_count;
		m_bound.resize(bound_count);
		for (index_type i = 0; i < bound_count; i++)
		{
			file >> m_bound[i].idx >> m_bound[i].node_type >> m_bound[i].ref_node;
			file >> m_bound[i].bound_normal[0] >> m_bound[i].bound_normal[1] >> m_bound[i].bound_normal[2];
		}

	}

	void GridConvertPiflow2Flexible::ReadCell(std::ifstream& file)
	{
		index_type cell_count;
		file >> cell_count;
		m_cell.resize(cell_count);
		for (index_type i = 0; i < cell_count; i++)
		{
			for (index_type j = 0; j < 8; j++)
			{
				file >> m_cell[i].node[j];
			}
		}

	}

	void GridConvertPiflow2Flexible::WriteNodeCoord(std::ofstream& file)
	{
		file << m_coord.size() << std::endl;
		for (index_type i = 0; i < m_coord.size(); i++)
		{
			file << m_coord[i].x << " " << m_coord[i].y << " " << m_coord[i].z << std::endl;
		}

	}

	void GridConvertPiflow2Flexible::WriteInnerNode(std::ofstream& file)
	{
		file << m_inner_node.size() << std::endl;
		for (index_type i = 0; i < m_inner_node.size(); i++)
		{
			file << m_inner_node[i].idx << " ";
			for (index_type j = 0; j < 6; j++)
			{
				file << m_inner_node[i].neighbor[j] << " ";
			}
			file << std::endl;
		}

	}

	void GridConvertPiflow2Flexible::WriteBound(std::ofstream& file)
	{
		file << m_bound.size() << std::endl;
		for (index_type i = 0; i < m_bound.size(); i++)
		{
			file << m_bound[i].idx << " " << m_bound[i].node_type << " " << m_bound[i].ref_node << " ";
			file << m_bound[i].bound_normal[0] << " " << m_bound[i].bound_normal[1] << " " << m_bound[i].bound_normal[2] << std::endl;
		}

	}

	void GridConvertPiflow2Flexible::WriteCell(std::ofstream& file)
	{
		file << m_cell.size() << std::endl;
		for (index_type i = 0; i < m_cell.size(); i++)
		{
			for (index_type j = 0; j < 8; j++)
			{
				file << m_cell[i].node[j] << " ";
			}
			file << std::endl;
		}
	}

}