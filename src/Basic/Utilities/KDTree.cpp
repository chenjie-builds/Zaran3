/*
 * file: KDTree.hpp
 * author: J. Frederico Carvalho
 *
 * This is an adaptation of the KD-tree implementation in rosetta code
 * https://rosettacode.org/wiki/K-d_tree
 *
 * It is a reimplementation of the C code using C++.  It also includes a few
 * more queries than the original, namely finding all points at a distance
 * smaller than some given distance to a point.
 *
 */

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <limits>
#include <memory>
#include <vector>

#include "KDTree.h"

KDNode::KDNode() = default;

KDNode::KDNode(const coord_vec & coord, const size_t & index, const KDNodePtr & leftNode, const KDNodePtr & rightNode)
{
	coord_ = coord;
	index_ = index;
	leftNode_ = leftNode;
	rightNode_ = rightNode;
}

KDNode::KDNode(const pointIndex & pi, const KDNodePtr & leftNode, const KDNodePtr & rightNode)
{
	coord_ = pi.first;
	index_ = pi.second;
	leftNode_ = leftNode;
	rightNode_ = rightNode;
}

KDNode::~KDNode() = default;

double KDNode::GetCoord(const size_t & idx) { return coord_.at(idx); }
KDNode::operator bool() { return (!coord_.empty()); }
KDNode::operator coord_vec() { return coord_; }
KDNode::operator size_t() { return index_; }
KDNode::operator pointIndex() { return pointIndex(coord_, index_); }

KDNodePtr NewKDNodePtr()
{
	KDNodePtr mynode = std::make_shared< KDNode >();
	return mynode;
}

inline double Dist2(const coord_vec & a, const coord_vec & b)
{
	double distance = 0;
	for (size_t i = 0; i < a.size(); i++)
	{
		double di = a.at(i) - b.at(i);
		distance += di * di;
	}
	return distance;
}

inline double Dist2(const KDNodePtr & a, const KDNodePtr & b)
{
	return Dist2(a->coord_, b->coord_);
}

inline double Dist(const coord_vec & a, const coord_vec & b)
{
	return std::sqrt(Dist2(a, b));
}

inline double Dist(const KDNodePtr & a, const KDNodePtr & b)
{
	return std::sqrt(Dist2(a, b));
}

Comparer::Comparer(size_t index) : idx_{ index } {};

inline bool Comparer::CompareIdx(const pointIndex & a, const pointIndex & b)
{
	return (a.first.at(idx_) < b.first.at(idx_));
}

inline void SortOnIdx(const pointIndexArr::iterator & begin, const pointIndexArr::iterator & end, size_t idx)
{
	Comparer comp(idx);
	comp.idx_ = idx;

	using std::placeholders::_1;
	using std::placeholders::_2;

	std::nth_element(begin, begin + std::distance(begin, end) / 2, end, std::bind(&Comparer::CompareIdx, comp, _1, _2));
}

using point_vec = std::vector< coord_vec >;

KDNodePtr KDTree::MakeTree(const pointIndexArr::iterator & begin, const pointIndexArr::iterator & end, const size_t & length, const size_t & level)
{
	if (begin == end)
	{
		return NewKDNodePtr();  // empty tree
	}
	size_t dim = begin->first.size();
	
	//¸ù¾ÝlevelÅÅÐò
	if (length > 1)
	{
		SortOnIdx(begin, end, level);
	}

	auto middle = begin + (length / 2);

	auto l_begin = begin;
	auto l_end = middle;
	auto r_begin = middle + 1;
	auto r_end = end;

	size_t l_len = length / 2;
	size_t r_len = length - l_len - 1;

	KDNodePtr left;
	if (l_len > 0 && dim > 0)
	{
		left = MakeTree(l_begin, l_end, l_len, (level + 1) % dim);
	}
	else
	{
		left = leaf_;
	}
	KDNodePtr right;
	if (r_len > 0 && dim > 0)
	{
		right = MakeTree(r_begin, r_end, r_len, (level + 1) % dim);
	}
	else
	{
		right = leaf_;
	}

	// KDNode result = KDNode();
	return std::make_shared< KDNode >(*middle, left, right);
}

KDTree::KDTree(point_vec point_array)
{
	leaf_ = std::make_shared< KDNode >();
	// iterators
	pointIndexArr arr;
	for (size_t i = 0; i < point_array.size(); i++)
	{
		arr.push_back(pointIndex(point_array.at(i), i));
	}

	auto begin = arr.begin();
	auto end = arr.end();

	size_t length = arr.size();
	size_t level = 0;  // starting

	root_ = KDTree::MakeTree(begin, end, length, level);
}

KDTree::KDTree(point_vec point_array, index_vec index_array)
{

	leaf_ = std::make_shared< KDNode >();
	// iterators
	pointIndexArr arr;
	for (size_t i = 0; i < point_array.size(); i++)
	{
		arr.push_back(pointIndex(point_array.at(i), index_array.at(i)));
	}

	auto begin = arr.begin();
	auto end = arr.end();

	size_t length = arr.size();
	size_t level = 0;  // starting

	root_ = KDTree::MakeTree(begin, end, length, level);
}

KDTree::KDTree(pointIndexArr point_index_array)
{
	leaf_ = std::make_shared< KDNode >();
	auto begin = point_index_array.begin();
	auto end = point_index_array.end();

	size_t length = point_index_array.size();
	size_t level = 0;  // starting

	root_ = KDTree::MakeTree(begin, end, length, level);

}

KDNodePtr KDTree::Nearest(const KDNodePtr & branch, const coord_vec & pt, const size_t & level, const KDNodePtr & best, const double& best_dist)const
{
	double d, dx, dx2;

	if (!bool(*branch))
	{
		return NewKDNodePtr();  // basically, null
	}

	coord_vec branch_pt(*branch);
	size_t dim = branch_pt.size();

	d = Dist2(branch_pt, pt);
	dx = branch_pt.at(level) - pt.at(level);
	dx2 = dx * dx;

	KDNodePtr best_l = best;
	double best_dist_l = best_dist;

	if (d < best_dist)
	{
		best_dist_l = d;
		best_l = branch;
	}

	size_t next_lv = (level + 1) % dim;
	KDNodePtr section;
	KDNodePtr other;

	// select which branch makes sense to check
	if (dx > 0)
	{
		section = branch->leftNode_;
		other = branch->rightNode_;
	}
	else
	{
		section = branch->rightNode_;
		other = branch->leftNode_;
	}

	// keep nearest neighbor from further down the tree
	KDNodePtr further = Nearest(section, pt, next_lv, best_l, best_dist_l);
	if (!further->coord_.empty())
	{
		double dl = Dist2(further->coord_, pt);
		if (dl < best_dist_l)
		{
			best_dist_l = dl;
			best_l = further;
		}
	}
	// only check the other branch if it makes sense to do so
	if (dx2 < best_dist_l)
	{
		further = Nearest(other, pt, next_lv, best_l, best_dist_l);
		if (!further->coord_.empty())
		{
			double dl = Dist2(further->coord_, pt);
			if (dl < best_dist_l)
			{
				best_dist_l = dl;
				best_l = further;
			}
		}
	}

	return best_l;
};

// default caller
KDNodePtr KDTree::Nearest(const coord_vec & pt)const
{
	size_t level = 0;
	// KDNodePtr best = branch;
	double branch_dist = Dist2(coord_vec(*root_), pt);
	return Nearest(root_,          // beginning of tree
		pt,            // point we are querying
		level,         // start from level 0
		root_,          // best is the root
		branch_dist);  // best_dist = branch_dist
};

coord_vec KDTree::NearestPoint(const coord_vec & pt)const
{
	return coord_vec(*Nearest(pt));
};
size_t KDTree::NearestIndex(const coord_vec & pt)const
{
	return size_t(*Nearest(pt));
};

pointIndex KDTree::NearestPointIndex(const coord_vec & pt)const
{
	KDNodePtr nearest = Nearest(pt);
	return pointIndex(coord_vec(*nearest), size_t(*nearest));
}

pointIndexArr KDTree::Neighborhood(const KDNodePtr & branch, const coord_vec & pt, const double& rad, const size_t & level)
{
	double d, dx, dx2;

	if (!bool(*branch))
	{
		// branch has no point, means it is a leaf,
		// no points to add
		return pointIndexArr();
	}

	size_t dim = pt.size();

	double r2 = rad * rad;

	d = Dist2(coord_vec(*branch), pt);
	dx = coord_vec(*branch).at(level) - pt.at(level);
	dx2 = dx * dx;

	pointIndexArr nbh, nbh_s, nbh_o;
	if (d <= r2)
	{
		nbh.push_back(pointIndex(*branch));
	}

	//
	KDNodePtr section;
	KDNodePtr other;
	if (dx > 0)
	{
		section = branch->leftNode_;
		other = branch->rightNode_;
	}
	else
	{
		section = branch->rightNode_;
		other = branch->leftNode_;
	}

	nbh_s = Neighborhood(section, pt, rad, (level + 1) % dim);
	nbh.insert(nbh.end(), nbh_s.begin(), nbh_s.end());
	if (dx2 < r2)
	{
		nbh_o = Neighborhood(other, pt, rad, (level + 1) % dim);
		nbh.insert(nbh.end(), nbh_o.begin(), nbh_o.end());
	}

	return nbh;
};

pointIndexArr KDTree::Neighborhood(const coord_vec & pt, const double& rad)
{
	size_t level = 0;
	return Neighborhood(root_, pt, rad, level);
}

point_vec KDTree::NeighborhoodPoints(const coord_vec & pt, const double& rad)
{
	size_t level = 0;
	pointIndexArr nbh = Neighborhood(root_, pt, rad, level);
	point_vec nbhp;
	nbhp.resize(nbh.size());
	std::transform(nbh.begin(), nbh.end(), nbhp.begin(), [](pointIndex x) { return x.first; });
	return nbhp;
}

index_vec KDTree::NeighborhoodIndices(const coord_vec & pt, const double& rad)
{
	size_t level = 0;
	pointIndexArr nbh = Neighborhood(root_, pt, rad, level);
	index_vec nbhi;
	nbhi.resize(nbh.size());
	std::transform(nbh.begin(), nbh.end(), nbhi.begin(), [](pointIndex x) { return x.second; });
	return nbhi;
}
