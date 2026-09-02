/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DEMDataManager.h
 * \brief DEMDataManager: manages DEMFieldData (no grid-based FieldData needed).
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include "BasicType.h"
#include "DEMFieldData.h"
namespace zaran
{
    /// @brief DEM 数据管理器，直接持有 DEMFieldData（不使用基于节点/单元的 FieldData）
    class DEMDataManager
    {
    public:
        explicit DEMDataManager(shared_ptr<DEMFieldData> dem_data)
            : m_dem_data(std::move(dem_data)) {}
        ~DEMDataManager() = default;

        shared_ptr<DEMFieldData> GetDEMData() const { return m_dem_data; }

    private:
        shared_ptr<DEMFieldData> m_dem_data;
    };
} // namespace zaran
