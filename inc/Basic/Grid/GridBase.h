/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file GridBase.h
 * \brief Base class for grid.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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
		Block,
		Unkown,
    };
	class GridBase
	{
	public:
		GridBase(const string& name, index_type index, dimension_type dim,GridType type)
			:m_name(name),m_index(index),m_dim(dim),m_type(type){}
		virtual ~GridBase() {}
		void SetName(const std::string& name) { m_name = name; }
		void SetIndex(const index_type& index) { m_index = index; }
		void SetDimension(const index_type& dim) { m_dim = dim; }
		void SetType(GridType type = GridType::Unkown) { m_type = type; };
	public:
		const string& GetName()const { return m_name; }
		const index_type& GetIndex()const { return m_index; }
		const dimension_type& GetDim()const { return m_dim; }
		const GridType& GetType()const { return m_type; }
	private:
		string m_name;//网格名称
		index_type m_index;//网格索引
		dimension_type m_dim;//网格维度
		GridType m_type;//网格类型
	};
}