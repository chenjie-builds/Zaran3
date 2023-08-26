#include"StarCDGrid.h"
#include<fstream>
#include<iostream>
StarCDGrid::StarCDGrid()
{
	dim_ = Dimension::two;
	nodeNum_ = 0;
	elementNum_ = 0;
	faceNum_ = 0;
	triaNum_ = quadNum_ = tetrNum_ = 0;
	pyraNum_ = prsmNum_ = hexaNum_ = 0;
	face2Num_ = face3Num_ = face4Num_ = 0;
	element_.resize(0);
	boundFace_.resize(0);
	coord_.resize(0);//节点坐标数组

}
void StarCDGrid::ReadGridFromStraCD(string& filename)
{
	read_inp(filename + ".inp");
	read_vrt(filename + ".vrt");
	read_cel(filename + ".cel");
	read_bnd(filename + ".bnd");

}

vector<int> StarCDGrid::GetBoundFaceType() const
{
	vector <int> bound_type(boundFace_.size());
	for (size_t iFace = 0; iFace < boundFace_.size(); ++iFace)
	{
		bound_type[iFace] = boundFace_[iFace].GetBoundType();
	}
	return bound_type;
}
int StarCDGrid::GetBoundFaceNum() const
{
	return boundFace_.size();
}
vector<vector<int>> StarCDGrid::GetNodeInBoundFaceIndex() const
{
	vector<vector<int>> nodeVec(boundFace_.size());
	for (size_t iFace = 0; iFace < boundFace_.size(); ++iFace)
	{
		nodeVec[iFace] = boundFace_[iFace].GetNodeIndexVec();
	}
	return nodeVec;

}
vector<vector<int>> StarCDGrid::GetNodeInElementIndex() const
{

	vector<vector<int>> nodeVec(element_.size());
	for (size_t iElem = 0; iElem < element_.size(); ++iElem)
	{
		nodeVec[iElem] = element_[iElem].GetNodeIndexVec();
	}
	return nodeVec;

}
void StarCDGrid::read_inp(string& filename)
{
	std::ifstream fin(filename);
	string line;
	while (!fin.eof())
	{
		std::getline(fin, line);
		if (line.substr(0, 4) == "CTAB" || line.substr(0, 4) == "ctab")
		{

		}
		else if (line.substr(0, 5) == "RDEF" || line.substr(0, 5) == "rdef")
		{

		}
		else if (line.substr(0, 4) == "RNAME" || line.substr(0, 4) == "rname")
		{

		}
	}
	fin.close();
}

void StarCDGrid::read_vrt(string& filename)
{
	std::ifstream fin(filename);
	double ip, x, y, z;
	while (fin >> ip || !fin.eof())
	{
		fin >> x;
		fin >> y;
		fin >> z;
		coord_.emplace_back(Eigen::Vector3d(x, y, z));
	}
	nodeNum_ = coord_.size();
	fin.close();
}

void StarCDGrid::read_cel(string& filename)
{
	vector<int>point;
	int blockIndex;
	int id;
	std::ifstream fin(filename);
	if (dim_ == Dimension::two)
	{
		point.resize(6);
		while (fin >> id || !fin.eof())
		{

			fin >> point[0];
			fin >> point[1];
			fin >> point[2];
			fin >> point[3];
			fin >> blockIndex;

			fin >> point[5];
			if (point[2] == point[3])
			{
				element_.push_back(Element(vector<int>{point[0] - 1, point[1] - 1, point[2] - 1}, blockIndex));
				++triaNum_;
			}
			else
			{
				element_.push_back(Element(vector<int>{point[0] - 1, point[1] - 1, point[2] - 1, point[3] - 1}, blockIndex));
				quadNum_++;
			}
			if (blockIndex_.size() == 0)
			{
				blockIndex_.emplace_back(blockIndex);
				continue;
			}
			for (size_t iBlock = 0; iBlock < blockIndex_.size(); ++iBlock)
			{
				if (blockIndex == blockIndex_[iBlock])
					break;
				else if (iBlock == blockIndex_.size() - 1)
				{
					blockIndex_.emplace_back(blockIndex);
					break;
				}
			}
		}
		elementNum_ = triaNum_ + quadNum_;
	}
	fin.close();
}

void StarCDGrid::read_bnd(string& filename)
{
	std::ifstream fin(filename);
	vector<int>point;
	int id, iFace, bound_type;

	if (dim_ == Dimension::two)
	{
		point.resize(4);
		while (fin >> id || !fin.eof())
		{
			fin >> point[0];
			fin >> point[1];
			fin >> point[2];
			fin >> point[3];
			fin >> iFace;
			fin >> iFace;
			fin >> bound_type;
			if (point[2] == point[3])
			{
				face2Num_++;
				boundFace_.push_back(BoundFace(vector<int>{point[0] - 1, point[1] - 1}, bound_type));
			}

		}
	}
}
StarCDGrid::Element::Element()
{
	nodeIndexVec_.resize(0);
	blockIndex_ = -1;
}
StarCDGrid::Element::Element(const vector<int>& node_index, const int& block_index)
{
	SetNodeIndex(node_index);
	SetBlockIndex(block_index);
}
StarCDGrid::BoundFace::BoundFace()
{
	nodeIndexVec_.resize(0);
	boundType_ = -1;
}
StarCDGrid::BoundFace::BoundFace(const vector<int>& node_index, const int& bound_type)
{
	SetNodeIndex(node_index);
	SetType(bound_type);
}