//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	StarCDGrid.h														||
//*	@brief	StarCD 网格类, 用于读取starCD 文件格式								||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"CommonPara.h"
#include <Eigen/Dense>
#include<string>
#include<vector>
using std::string;
using std::vector;
class StarCDGrid
{
public:
	StarCDGrid();
	void ReadGridFromStraCD(string& filename);
public:
	// get the dimension of grid
	Dimension GetDimension()const { return dim_; }
	// get the point number
	int GetNodeNum()const { return nodeNum_; }
	// get the element number
	int GetElementNum()const { return elementNum_; }
	// get the face number
	int GetFaceNum()const { return faceNum_; }

	// get the triangle element number
	int GetTriaNum()const { return triaNum_; }
	// get the quadrangle element number
	int GetQuadNum()const { return quadNum_; }
	// get the tetrahedron element number
	int GetTetrNum()const { return tetrNum_; }
	// get the pyramid element number
	int GetPyraNum()const { return pyraNum_; }
	// get the prism element number
	int GetPrsmNum()const { return prsmNum_; }
	// get the hexahedron element number
	int GetHexaNum()const { return hexaNum_; }

	// get the face number with two nodes
	int GetFaceNum2()const { return face2Num_; }
	// get the face number with three nodes
	int GetFaceNum3()const { return face3Num_; }
	// get the face number with four nodes
	int GetFaceNum4()const { return face4Num_; }

	// get the nodes' coordinate
	vector<Eigen::Vector3d>& GetNodeCoord() { return coord_; }
	// get the boundary face type
	vector<int> GetBoundFaceType()const;
	// get the boundary face number
	int GetBoundFaceNum()const;
	// get the nodes' index in boundary faces
	vector<vector<int>> GetNodeInBoundFaceIndex()const;
	// get the nodes' index in elements
	vector<vector<int>> GetNodeInElementIndex()const;
	// get the block number: the last element's block index
	int GetBlockNum()const { return blockIndex_.size(); }

	vector<int> GetBlockIndex()const { return blockIndex_; }
	// gte the block index
	int GetElementBlock(size_t iElem)const { return element_[iElem].GetBolockIndex(); }
private:
	class Element
	{
	public:
		Element();
		Element(const vector<int>& node_index, const int& blockIndex);
	public:
		void SetNodeIndex(const vector<int>& node_index) { nodeIndexVec_ = node_index; }
		void SetBlockIndex(const int& block_index) { blockIndex_ = block_index; }
	public:
		vector<int>GetNodeIndexVec() const { return nodeIndexVec_; }
		int GetNodeIndexVec(size_t iNode)const { return nodeIndexVec_[iNode]; }
		int GetNodeNum() const { return nodeIndexVec_.size(); }
		int GetBolockIndex() const { return blockIndex_; }
	private:
		vector<int>nodeIndexVec_;
		int blockIndex_;
	};
	class BoundFace
	{
	public:
		BoundFace();
		BoundFace(const vector<int>& node_index, const int& bound_type);
	public:
		void SetNodeIndex(const vector<int>& node_index) { nodeIndexVec_ = node_index; }
		void SetType(const int& bound_type) { boundType_ = bound_type; }
	public:
		vector<int>GetNodeIndexVec() const { return nodeIndexVec_; }
		int GetNodeIndexVec(size_t iNode)const { return nodeIndexVec_[iNode]; }
		int GetNodeNum() const { return nodeIndexVec_.size(); }
		int GetBoundType() const { return boundType_; }
	private:
		vector<int>nodeIndexVec_;
		int boundType_;
	};
private:
	void read_inp(string& filename);
	void read_vrt(string& filename);
	void read_cel(string& filename);
	void read_bnd(string& filename);

private:
	Dimension dim_;
	int nodeNum_, elementNum_, faceNum_;
	int triaNum_, quadNum_, tetrNum_;
	int pyraNum_, prsmNum_, hexaNum_;
	int face2Num_, face3Num_, face4Num_;

	vector<Element> element_;
	vector<BoundFace> boundFace_;
	vector<Eigen::Vector3d> coord_;//节点坐标数组
	vector<int> blockIndex_;
};