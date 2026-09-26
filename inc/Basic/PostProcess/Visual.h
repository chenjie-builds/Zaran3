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

		/// @brief 将 DEM 粒子与键合输出为 Tecplot ASCII (.dat)
		/// @details Tecplot ASCII 的文件头 `VARIABLES` 是**整个文件共享**的（不是每个 zone
		///          各自一套）。粒子是点云（约 25~45 个变量），键合是线元（3 个节点量 +
		///          约 20 个单元中心量），两者变量数不同，写进同一文件只能用大量空列去凑，
		///          因此拆成两个文件：
		///            - `result/particles_<iter>.dat`：ZONETYPE=FEPOINT、DATAPACKING=POINT，
		///              一行一个粒子，便于脚本直接解析；
		///            - `result/bonds_<iter>.dat`：ZONETYPE=FELINESEG、DATAPACKING=BLOCK，
		///              键合量以 VARLOCATION=CELLCENTERED 输出（Tecplot 要求带单元中心量的
		///              zone 必须用 BLOCK，且节点变量块在前、单元中心变量块在后、连接表最后）。
		/// @param field_manager 场管理器（自动寻找第一个 DEMField）
		/// @param iter 当前迭代步（用于文件命名）
		/// @param solution_time 当前物理时间，写入 ZONE 的 SOLUTIONTIME 供 Tecplot 做时间动画
		/// @param with_bonds 是否同时写键合文件（对应 output.bond_details；键合数为 0 时自动跳过）
		void WriteParticleTecplotASCII(shared_ptr<FieldManager> field_manager, int iter,
			double solution_time, bool with_bonds);

		/// @brief 数据级接口：直接由 DEMFieldData 写两个 Tecplot ASCII 文件
		/// @param particle_file 粒子文件全路径
		/// @param bond_file 键合文件全路径（with_bonds 为 false 或没有键合时不创建）
		void WriteParticleTecplotASCII(const DEMFieldData& dem_data,
			const std::string& particle_file, const std::string& bond_file,
			double solution_time, bool with_bonds);

		/// @brief 均匀笛卡尔网格上的一个标量场（供"非求解器"算例如粒子→网格映射使用）
		struct UniformScalarField
		{
			std::string name;                              ///< Tecplot 变量名
			const std::vector<double>* values = nullptr;   ///< 长度 = ni*nj*nk，k 最慢、i 最快
		};

		/// @brief 把均匀笛卡尔网格上的若干标量场写成 Tecplot ASCII（ORDERED + POINT）
		/// @details 与网格文件无关：网格由"首节点坐标 + 间距 + 节点数"直接给出，
		///          与 ParticleGridMapper 的节点型控制体约定一致
		///          （首节点是第 0 个单元的中心，计算域两端各多出半格）。
		///          变量表恒为 X,Y,Z,<fields...>。
		/// @param zone_name 写进 ZONE 的 T=，建议带上帧的说明（如 "alpha subcell n=8"），
		///                  Tecplot 里一眼就能分辨是哪一帧
		void WriteUniformGridScalarTecplotASCII(const std::string& filename,
			const std::string& title, const std::string& zone_name,
			int ni, int nj, int nk,
			double x0, double y0, double z0,
			double dx, double dy, double dz,
			const std::vector<UniformScalarField>& fields,
			double solution_time);

		/// @brief 点云上的一个标量场（供相间力等"自定义列"输出使用）
		struct PointScalarField
		{
			std::string name;                              ///< Tecplot 变量名
			const std::vector<double>* values = nullptr;   ///< 长度 = 点数
		};

		/// @brief 把一批点（及其上的若干标量）写成 Tecplot ASCII 点云
		/// @details `ZONETYPE=FEPOINT` + `DATAPACKING=POINT`，变量表 = 各场的名字，
		///          一行一个点，脚本可直接按行解析（与 DEM 的 particles_*.dat 同口径）。
		///          每 512 个值强制换行：Tecplot ASCII 单行上限 32000 字节，超了会静默截断。
		///          DEM 的 `WriteParticleTecplotASCII` 列是固定的，这里是"任意列"的通用版本。
		void WritePointsTecplotASCII(const std::string& filename,
			const std::string& title, const std::string& zone_name,
			int num_points, const std::vector<PointScalarField>& fields,
			double solution_time);
	};
}