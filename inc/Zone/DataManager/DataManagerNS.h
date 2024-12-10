#pragma once
#include"FieldData.h"
#include "GridFNFDM.h"
#include "FieldDataManager.h"
namespace zaran
{
    /// @brief NS方程数据管理类
    /// @details 用于管理NS方程的数据
    class DataManagerNS :public DataManager
    {
    public:
        DataManagerNS(std::shared_ptr<FieldData> fieldData, int data_num);
        ~DataManagerNS();
    public:
        void CreateData()override;
        void RegisterData()override;
        int GetEqNum() const { return m_equ_num; }
    public:
        //! @brief 设置基本变量
        /// @param idx_eq 变量索引
        /// @param idx_data 数据索引 
        /// @param value 输入的值 
        void SetPrim(int idx_eq, int idx_data, double value);
		void SetPrim(int idx_data, const double* value);
        //! @brief 设置守恒变量
        /// @param idx_eq 变量索引
        /// @param idx_data 数据索引
        /// @param value 输入的值
        void SetCons(int idx_eq, int idx_data, double value);
		void SetCons(int idx_data, const double* value);
        //! @brief 设置上一个时间步守恒变量
        /// @param idx_eq 变量索引
        /// @param idx_data 数据索引
        /// @param value 输入的值
        void SetConsOld(int idx_eq, int idx_data, double value);
        //! @brief 设置残差
        /// @param idx_eq 变量索引
        /// @param idx_data 数据索引
        /// @param value 输入的值
        void SetResidual(int idx_eq, int idx_data, double value);
        void SetResidual(int idx_data, const double* value);
        //! @brief 设置限制器系数
        /// @param idx_eq 变量索引
        /// @param idx_data 数据索引
        /// @param value 输入的值
        void SetLimiter(int idx_eq, int idx_data, double value);
        //! @brief 设置粘性通量
        /// @param idx_eq 变量索引
        /// @param idx_dim 通量维度索引
        /// @param idx_data 数据索引
        /// @param value 输入的值
        void SetViscousFlux(int idx_eq, int idx_dim, int idx_data, double value);
        //! @brief 设置粘性通量梯度
        /// @param idx_eq 变量索引
        /// @param idx_flux 通量索引, 0:x方向, 1:y方向, 2:z方向
        /// @param idx_grad_dim 梯度维度索引
        /// @param idx_data 数据索引
        /// @param value 输入的值
        void SetViscousFluxGrad(int idx_eq, int idx_flux, int idx_grad_dim, int idx_data, double value);
        //! @brief 设置原始变量梯度
        /// @param idx_eq 变量索引
        /// @param idx_grad_dim 梯度维度索引
        /// @param idx_data 数据索引
        /// @param value 输入的值
        void SetPrimitiveGrad(int idx_eq, int idx_grad_dim, int idx_data, double value);
        //! @brief 设置温度
        /// @param idx_data 数据索引
        /// @param value 输入的值
        void SetTemperature(int idx_data, double value);
        //! @brief 设置温度梯度
        /// @param idx_grad_dim 梯度维度索引
        /// @param idx_data 数据索引
        void SetTemperatureGrad(int idx_grad_dim, int idx_data, double value);
        //! @brief 设置非物理点标记
        /// @param idx_data 数据索引
        /// @param value 输入的值
        void SetNonPhysical(int idx_data, int value);
        //! @brief 设置时间步长
        /// @param idx_data 数据索引
        /// @param dt 输入的值
        void SetTimeStep(int idx_data, double dt);
    public:
        //! @brief 获取基本变量
        /// @param idx_eq 变量索引
        /// @param idx_data 数据索引
        double GetPrim(int idx_eq, int idx_data);
        //! @brief 获取基本变量
        /// @param idx_eq 变量索引
        /// @return 第 idx_eq 个基本变量
        double* GetPrim(int idx_eq);
        double* GetDensity();
        double GetDensity(int idx_data);
        double* GetVelocity(int idx_dim);
        //! @brief 获取速度
        /// @param idx_dim 速度维度索引
        /// @param idx_data 数据索引
        double GetVelocity(int idx_dim, int idx_data);
        double* GetPressure();
        double GetPressure(int idx_data);
        double GetTimeStep(int idx_data);
        double* GetTemperature();
        double GetTemperature(int idx_data);
        double* GetConservative(int idx_eq);
        double* GetConsOld(int idx_eq);
        double GetCons(int idx_eq, int idx_data);
        double GetConsOld(int idx_eq, int idx_data);
        double* GetResidual(int idx_eq);
        double GetResidual(int idx_eq, int idx_data);
        double* GetLimiter(int idx_eq);
        double GetLimiter(int idx_eq, int idx_data);
        double* GetViscousFlux(int idx_eq, int idx_dim);
        double GetViscousFlux(int idx_eq, int idx_dim, int idx_data);
        double* GetViscousFluxGrad(int idx_eq, int idx_dim, int idx_grad_dim);
        double GetViscousFluxGrad(int idx_eq, int idx_dim, int idx_grad_dim, int idx_data);
        double* GetPrimitiveGrad(int idx_eq, int idx_grad_dim);
        double GetPrimGrad(int idx_eq, int idx_grad_dim, int idx_data);
        int* GetNonPhysical();
        int GetNonPhysical(int idx_data);
    protected:
        int m_equ_num = 5;
        /// @brief 基本变量, 1D 变量索引, 2D 数据索引
        double** m_prim;
        /// @brief 守恒变量, 1D 变量索引, 2D 数据索引
        double** m_cons;
        /// @brief 上一个时间步原始变量, 用于Runge-Kutta, 1D 变量索引, 2D 数据索引
        double** m_cons_old;
        /// @brief 时间步长
        double* m_dt;
        /// @brief 残差, 1D 变量索引, 2D 数据索引
        double** m_residual;
        /// @brief 限制器, 1D 变量索引, 2D 数据索引
        double** m_limiter;
        /// @brief 温度
        double* m_temperture;
        /// @brief 温度梯度, 1D 梯度维度索引, 2D 数据索引
        double** m_temperture_grad;
        /// @brief 基本变量梯度, 1D 变量索引, 2D 梯度维度索引, 3D 数据索引
        double*** m_prim_grad;
        /// @brief 粘性通量, 1D 变量索引, 2D 通量维度索引, 3D 数据索引
        double*** m_viscous_flux;
        /// @brief 粘性通量梯度, 1D 变量索引, 2D 通量维度索引, 3D 梯度维度索引, 4D 数据索引
        double**** m_viscous_flux_grad;
        /// @brief 非物理点标记
        int* m_non_physical;
    };

}