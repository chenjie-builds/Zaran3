//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	GridBase.h															||
//*	@brief	网格基类																||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include"BasicType.h"
#include "CommonPara.h"
namespace zaran
{
    enum class GridType
    {
        Structured,
        Flexible,
        Zaran,
		Unkown,
    };
	class GridBase
	{
	public:
		GridBase() :name_(), index_(0), dim_(), type_(GridType::Unkown) {}
		virtual ~GridBase() {}
		void SetName(const std::string& name) { name_ = name; }
		void SetIndex(const size_t& index) { index_ = index; }
		void SetDimension(const Dimension& dim) { dim_ = dim; }
		void SetType(GridType type = GridType::Unkown) { type_ = type; };
	public:
		const string& GetName()const { return name_; }
		const size_t& GetIndex()const { return index_; }
		const Dimension& GetDimension()const { return dim_; }
		const GridType& GetType()const { return type_; }
	private:
		string name_;//网格名称
		size_t index_;//网格索引
		Dimension dim_;//网格维度
		GridType type_;//网格类型
	};
}