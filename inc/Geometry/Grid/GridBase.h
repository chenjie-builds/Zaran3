//==============================================================================||
//*	ZaRan	-	A Totally Automatic CFD Software								||
//*	Copyright (C) ,Since 2020													||
//*-----------------------------------------------------------------------------||
//*	License																		||
//*	This file is part of ZaRan.													||
//*																				||
//*	@file	GridBase.h															||
//*	@brief	Íø¸ñ»ùÀà																||
//*	@author	Chen Jie.															||
//==============================================================================||
#pragma once
#include <string>
#include "CommonPara.h"
namespace zaran
{
	enum class GridType
	{
		STRUCT = 0,
		UNSTRUCT = 1,
		NOTSET = 100
	};
	class GridBase
	{
	public:
		GridBase() :name_(), index_(0), dim_(), type_(GridType::NOTSET) {}
		void SetName(const std::string& name) { name_ = name; }
		void SetIndex(const size_t& index) { index_ = index; }
		void SetDimension(const Dimension& dim) { dim_ = dim; }
		void SetType(GridType type = GridType::NOTSET) { type_ = type; };
	public:
		const std::string& GetName()const { return name_; }
		const size_t& GetIndex()const { return index_; }
		const Dimension& GetDimension()const { return dim_; }
		const GridType& GetType()const { return type_; }
	private:
		std::string name_;
		size_t index_;
		Dimension dim_;
		GridType type_;
	};
}