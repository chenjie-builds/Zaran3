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
		GridBase(const string& name,int index,int dim,GridType type)
			:m_name(name),m_index(index),m_dim(dim),m_type(type){}
		virtual ~GridBase() {}
		void SetName(const std::string& name) { m_name = name; }
		void SetIndex(const size_t& index) { m_index = index; }
		void SetDimension(const int& dim) { m_dim = dim; }
		void SetType(GridType type = GridType::Unkown) { m_type = type; };
	public:
		const string& GetName()const { return m_name; }
		const size_t& GetIndex()const { return m_index; }
		const int& GetDimension()const { return m_dim; }
		const GridType& GetType()const { return m_type; }
	private:
		string m_name;//网格名称
		int m_index;//网格索引
		int m_dim;//网格维度
		GridType m_type;//网格类型
	};
}