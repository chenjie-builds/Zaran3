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
#include"BasicType.h"
#include"CommonPara.h"
namespace zaran
{

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
		Array<DVector3D>& GetNodeCoord() { return coord_; }
		// get the boundary face type
		IArray GetBoundFaceType()const;
		// get the boundary face number
		int GetBoundFaceNum()const;
		// get the nodes' index in boundary faces
		Array<IArray> GetNodeInBoundFaceIndex()const;
		// get the nodes' index in elements
		Array<IArray> GetNodeInElementIndex()const;
		// get the block number: the last element's block index
		int GetBlockNum()const { return blockIndex_.size(); }

		IArray GetBlockIndex()const { return blockIndex_; }
		// gte the block index
		int GetElementBlock(size_t iElem)const { return element_[iElem].GetBolockIndex(); }
	private:
		class Element
		{
		public:
			Element();
			Element(const IArray& node_index, const int& blockIndex);
		public:
			void SetNodeIndex(const IArray& node_index) { nodeIndexVec_ = node_index; }
			void SetBlockIndex(const int& block_index) { blockIndex_ = block_index; }
		public:
			IArray GetNodeIndexVec() const { return nodeIndexVec_; }
			int GetNodeIndexVec(size_t iNode)const { return nodeIndexVec_[iNode]; }
			int GetNodeNum() const { return nodeIndexVec_.size(); }
			int GetBolockIndex() const { return blockIndex_; }
		private:
			IArray nodeIndexVec_;
			int blockIndex_;
		};
		class BoundFace
		{
		public:
			BoundFace();
			BoundFace(const IArray& node_index, const int& bound_type);
		public:
			void SetNodeIndex(const IArray& node_index) { nodeIndexVec_ = node_index; }
			void SetType(const int& bound_type) { boundType_ = bound_type; }
		public:
			IArray GetNodeIndexVec() const { return nodeIndexVec_; }
			int GetNodeIndexVec(size_t iNode)const { return nodeIndexVec_[iNode]; }
			int GetNodeNum() const { return nodeIndexVec_.size(); }
			int GetBoundType() const { return boundType_; }
		private:
			IArray nodeIndexVec_;
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

		Array<Element> element_;
		Array<BoundFace> boundFace_;
		Array<DVector3D> coord_;//节点坐标数组
		IArray blockIndex_;
	};
}