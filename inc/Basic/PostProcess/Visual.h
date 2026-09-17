/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file Visual.h
 * \brief Visual class, used to output field data.
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include <cgnslib.h>
#include "Field.h"
#include "NSFieldStruct.h"
#include "NSFieldZaran.h"
#include "DEMFieldData.h"
namespace zaran
{
	class Visual
	{
	public:
		void WriteTecASCII(shared_ptr<FieldManager> field_manager);
		void WriteTecplotASCII(shared_ptr<NSFieldStruct> field, std::ostream& os);
		/// @brief 带气相体积分数的 Tecplot ASCII 输出（两相求解器使用）
		void WriteTecplotASCII(shared_ptr<NSFieldStruct> field, std::ostream& os,
			bool with_volume_fraction);
		static void WriteTecplotASCII(const shared_ptr<NSFieldZaran>& field, std::ostream& os);
        void WriteTecplotASCII(shared_ptr<NSFieldFNFDM> field, std::ostream &os);


		void WriteTecplotBinary(shared_ptr<FieldManager> field_manager);
		void WriteTecplotBinary(shared_ptr<NSFieldFNFDM> field);
		void WriteTecplotBinary(shared_ptr<NSFieldZaran> field);
		void WriteTecplotBinary(shared_ptr<NSFieldStruct> field);
		/// @brief 带气相体积分数的 Tecplot 二进制输出（两相求解器使用）
		void WriteTecplotBinary(shared_ptr<NSFieldStruct> field, bool with_volume_fraction);

		/// @brief 该场的数据里是否含气相体积分数（两相场才注册了这个量）
		static bool HasVolumeFraction(const shared_ptr<Field>& field);

		void WriteVtkASCII(shared_ptr<FieldManager> field_manager);
		void WriteVtkASCII(shared_ptr<NSFieldStruct> field, std::ostream& os);
		void WriteVtkASCII(shared_ptr<NSFieldZaran> field, std::ostream& os);
		void WriteVtkASCII(shared_ptr<NSFieldFNFDM> field, std::ostream& os);
		void WriteVtkBinary(shared_ptr<FieldManager> field_manager);
		void WriteVtkBinary(shared_ptr<NSFieldStruct> field, std::ostream& os);

		void WriteCGNS(shared_ptr<FieldManager> field_manager);
		void WriteCGNS(shared_ptr<NSFieldStruct> field, cgsize_t index_file, cgsize_t index_base);
		void WriteCGNS(shared_ptr<NSFieldZaran> field, cgsize_t index_file, cgsize_t index_base);
		void WriteCGNS(shared_ptr<NSFieldFNFDM> field, cgsize_t index_file, cgsize_t index_base);

		/// @brief 将 DEM 粒子数据输出为 VTK XML PolyData (.vtp) 文件
		/// @param field_manager 场管理器（会自动寻找 DEMField）
		/// @param iter 当前迭代步（用于文件命名）
		void WriteParticleVTP(shared_ptr<FieldManager> field_manager, int iter);

		/// @brief 将单个 DEMFieldData 输出为 .vtp 文件
		void WriteParticleVTP(const DEMFieldData& dem_data, const std::string& filename);
	};
}