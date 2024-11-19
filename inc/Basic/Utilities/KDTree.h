//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	KDtree.h															||
//*	@brief	kd树数据结构, 后续可能需要根据相关需求进行修改							||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once

/*
 * file: KDTree.h
 * author: J. Frederico Carvalho
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 *  https://rosettacode.org/wiki/K-d_tree
 * It is a reimplementation of the C code using C++.
 * It also includes a few more queries than the original
 *
 */

#include <algorithm>
#include <functional>
#include <memory>
#include <vector>

using coord_vec = std::vector< double >;
using index_vec = std::vector< size_t >;
using pointIndex = typename std::pair< std::vector< double >, size_t >;

class KDNode
{
public:
	using KDNodePtr = std::shared_ptr< KDNode >;
	size_t index_;
	coord_vec coord_;
	KDNodePtr leftNode_;
	KDNodePtr rightNode_;

	// initializer
	KDNode();
	KDNode(const coord_vec& coord, const size_t&index, const KDNodePtr&leftNode, const KDNodePtr&rightNode);
	KDNode(const pointIndex&, const KDNodePtr&, const KDNodePtr&);
	~KDNode();

	// getter
	double GetCoord(const size_t&);

	// conversions
	explicit operator bool();
	explicit operator coord_vec();
	explicit operator size_t();
	explicit operator pointIndex();
};

using KDNodePtr = std::shared_ptr< KDNode >;

KDNodePtr NewKDNodePtr();

// square euclidean distance
inline double Dist2(const coord_vec&, const coord_vec&);
inline double Dist2(const KDNodePtr&, const KDNodePtr&);

// euclidean distance
// 两个点之间的距离
inline double Dist(const coord_vec&, const coord_vec&);
inline double Dist(const KDNodePtr&, const KDNodePtr&);

// Need for sorting
class Comparer
{
public:
	size_t idx_;
	explicit Comparer(size_t idx);
	inline bool CompareIdx(const std::pair< std::vector< double >, size_t >&, const std::pair< std::vector< double >, size_t >&);
};

using pointIndexArr = typename std::vector< pointIndex >;

inline void SortOnIdx(const pointIndexArr::iterator&, const pointIndexArr::iterator&, size_t idx);

using point_vec = std::vector< coord_vec >;
class KDTree 
{
private:
	KDNodePtr root_;
	KDNodePtr leaf_;

	KDNodePtr MakeTree(const pointIndexArr::iterator& begin, const pointIndexArr::iterator& end, const size_t& length, const size_t& level);

public:
	KDTree() = default;
	explicit KDTree(point_vec point_array);
	explicit KDTree(point_vec point_array, index_vec index_array);
	explicit KDTree(pointIndexArr point_index_array);

private:
	KDNodePtr Nearest(const KDNodePtr& branch, const coord_vec& pt, const size_t& level, const KDNodePtr& best, const double& best_dist)const;

	// default caller
	KDNodePtr Nearest(const coord_vec& pt)const;

public:
	coord_vec NearestPoint(const coord_vec& pt)const;
	size_t NearestIndex(const coord_vec& pt)const;
	pointIndex NearestPointIndex(const coord_vec& pt)const;

private:
	pointIndexArr Neighborhood(const KDNodePtr& branch, const coord_vec& pt, const double& rad, const size_t& level);

public:
	pointIndexArr Neighborhood(const coord_vec& pt, const double& rad);

	point_vec NeighborhoodPoints(const coord_vec& pt, const double& rad);

	index_vec NeighborhoodIndices(const coord_vec& pt, const double& rad);
};