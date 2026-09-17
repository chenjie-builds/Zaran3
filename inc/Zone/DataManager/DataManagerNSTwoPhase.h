/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file DataManagerNSTwoPhase.h
 * \brief Data manager for the two-phase (volume-fraction) NS/Euler field.
 * \author Chen Jie.
 *
 *  在 DataManagerNSStruct 的基础上多注册一个节点量 "volume_fraction"（气相体积
 *  分数 ε_g）。走 FieldData 注册的用意有两个：
 *    1. 它是**流场数据**的一部分，应当随场一起被 Visual 输出到 Tecplot；
 *    2. 阶段 2 的粒子→网格体积分配只需往这个数组里写，不必改动求解器接口。
 */

#pragma once
#include "DataManagerNSStruct.h"

namespace zaran
{
	/// @brief 两相场的数据管理器：DataManagerNSStruct + 体积分数场。
	class DataManagerNSTwoPhase : public DataManagerNSStruct
	{
	public:
		DataManagerNSTwoPhase(shared_ptr<FieldData> field_data, int ni, int nj, int nk);
		~DataManagerNSTwoPhase();
	public:
		void CreateData() override;
		void RegisterData() override;
	public:
		/// @brief 气相体积分数场（节点中心，含 ghost）
		double* GetVolumeFraction() { return m_volume_fraction; }
		const double* GetVolumeFraction() const { return m_volume_fraction; }

	private:
		/// @brief 气相体积分数
		double* m_volume_fraction = nullptr;
	};
}
