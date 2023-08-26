#include"SYSU.h"
#include<fstream>
using namespace zaran;
zaran::GridListFactorySYSU::GridListFactorySYSU()
{
	m_fileName = "sysu.dat";
}
void zaran::GridListFactorySYSU::Create(std::shared_ptr<GridList>& gridList)
{
	ReadFile();
}

void zaran::GridListFactorySYSU::ReadFile()
{
	std::ifstream fin(m_fileName);
	fin >> m_NodeNum;
	m_NodeX.resize(m_NodeNum);
	m_NodeY.resize(m_NodeNum);
	m_NodeZ.resize(m_NodeNum);
	for (size_t i = 0; i < m_NodeNum; i++)
	{
		fin >> m_NodeX[i] >> m_NodeY[i] >> m_NodeZ[i];
	}
	m_BoundNodeNum = 0;
	double nBound;
	fin >> nBound;
	m_BoundNodeNum += nBound;
	int tempIndex1, tempIndex2;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> tempIndex1 >> tempIndex2;
		m_InletNodeIndex.push_back(tempIndex1);
		m_InletNeiborNodeIndex.push_back(tempIndex2);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> tempIndex1 >> tempIndex2;
		m_OutletNodeIndex.push_back(tempIndex1);
		m_OutletNeiborNodeIndex.push_back(tempIndex2);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> tempIndex1 >> tempIndex2;
		m_OutletNodeIndex.push_back(tempIndex1);
		m_OutletNeiborNodeIndex.push_back(tempIndex2);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> tempIndex1 >> tempIndex2;
		m_OutletNodeIndex.push_back(tempIndex1);
		m_OutletNeiborNodeIndex.push_back(tempIndex2);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> tempIndex1 >> tempIndex2;
		m_OutletNodeIndex.push_back(tempIndex1);
		m_OutletNeiborNodeIndex.push_back(tempIndex2);
	}
	fin >> nBound;
	m_BoundNodeNum += nBound;
	for (size_t i = 0; i < nBound; i++)
	{
		fin >> tempIndex1 >> tempIndex2;
		m_WallNodeIndex.push_back(tempIndex1);
		m_WallNeiborNodeIndex.push_back(tempIndex2);
	}

}
