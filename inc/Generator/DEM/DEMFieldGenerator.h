/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMFieldGenerator.h
 * \brief Factory for creating a DEMField from GlobalData parameters.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "GlobalField.h"
#include "DEMField.h"
namespace zaran
{
    /// @brief DEM 场生成器，根据控制文件创建并填充 DEMField
    class DEMFieldGenerator
    {
    public:
        DEMFieldGenerator() = default;
        ~DEMFieldGenerator() = default;

        /// @brief 创建 FieldManager（内含一个 DEMField）
        shared_ptr<FieldManager> Create();

    private:
        /// @brief 从 CSV 文件加载粒子
        void LoadParticlesFromFile(const shared_ptr<DEMField>& field) const;
        /// @brief 从包围盒范围随机生成粒子（CSV 文件不存在时使用）
        void GenerateParticlesFromBox(const shared_ptr<DEMField>& field) const;
        /// @brief 添加包围盒六面墙
        void AddBoxWalls(const shared_ptr<DEMField>& field) const;
    };
} // namespace zaran
