#include "ReadPLY.h"
#include "Log.h"
namespace zaran
{
	void PlyReader::ReadPlyFile(const char* filename)
	{
		Log::info("Start Read Ply file {}.", filename);
		std::ifstream file(filename);
		if (!file.is_open())
		{
			Log::error("Error: Cannot open file {}.", filename);
			return;
		}

		std::string line;
		bool header = true;
		int vertex_count = 0;

		// 解析头部
		while (std::getline(file, line))
		{
			if (line.find("end_header") != std::string::npos)
			{
				header = false;
				break;
			}
			if (header && line.find("element vertex") != std::string::npos)
			{
				vertex_count = std::stoi(line.substr(15));
				m_points.reserve(vertex_count);
			}
		}
		// 读取顶点数据
		while (std::getline(file, line))
		{
			// 如果已经读取了足够的顶点，则停止读取
			if (m_points.size() >= vertex_count)
			{
				break;
			}
			dynamic_array<double> p(3);
			std::istringstream iss(line);
			if (!(iss >> p[0] >> p[1] >> p[2]))
			{
				Log::error("Error: Invalid vertex data", filename);
				continue;
			}
			m_points.push_back(p);
		}
		file.close();
		Log::info("Read Ply file Succeed!");
		if (vertex_count != m_points.size())
		{
			Log::warn("Warning: Vertex count mismatch. Expected {}, got {}.", vertex_count, m_points.size());
		}
		Log::info("Points number: {}", vertex_count);
	}
}
