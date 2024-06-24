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
        DataManagerNS(FieldData* fieldData, int data_num);
        ~DataManagerNS();
    public:
        void CreateData()override;
        void RegisterData()override;
    public:
        //! @brief 设置基本变量
        /// @param iEqu 变量索引
        /// @param iNode 节点索引 
        /// @param value 输入的值 
        void SetPrim(int iEqu, int iNode, double value);
        //! @brief 设置守恒变量
        /// @param iEqu 变量索引
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetCons(int iEqu, int iNode, double value);
        //! @brief 设置残差
        /// @param iEqu 变量索引
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetResidual(int iEqu, int iNode, double value);
        //! @brief 设置限制器系数
        /// @param iEqu 变量索引
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetLimiter(int iEqu, int iNode, double value);
        //! @brief 设置粘性通量
        /// @param iEqu 变量索引
        /// @param iDim 通量维度索引
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetViscousFlux(int iEqu, int iDim, int iNode, double value);
        //! @brief 设置粘性通量梯度
        /// @param iEqu 变量索引
        /// @param iDim 通量维度索引
        /// @param iGradDim 梯度维度索引
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetViscousFluxGrad(int iEqu, int iDim, int iGradDim, int iNode, double value);
        //! @brief 设置原始变量梯度
        /// @param iEqu 变量索引
        /// @param iGradDim 梯度维度索引
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetPrimitiveGrad(int iEqu, int iGradDim, int iNode, double value);
        //! @brief 设置温度
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetTemperature(int iNode, double value);
        //! @brief 设置温度梯度
        /// @param iGradDim 梯度维度索引
        /// @param iNode 节点索引
        void SetTemperatureGrad(int iGradDim, int iNode, double value);
        //! @brief 设置非物理点标记
        /// @param iNode 节点索引
        /// @param value 输入的值
        void SetNonPhysical(int iNode, int value);
        //! @brief 设置时间步长
        /// @param iNode 节点索引
        /// @param dt 输入的值
        void SetTimeStep(int iNode, double dt);
    public:
        //! @brief 获取基本变量
        /// @param iEqu 变量索引
        /// @param iNode 节点索引
        double GetPrim(int iEqu, int iNode);
        //! @brief 获取基本变量
        /// @param iEqu 变量索引
        /// @return 第 iEqu 个基本变量
        double* GetPrim(int iEqu);
        double* GetDensity();
        double GetDensity(int iNode);
        double* GetVelocity(int iDim);
        //! @brief 获取速度
        /// @param iDim 速度维度索引
        /// @param iNode 节点索引
        double GetVelocity(int iDim, int iNode);
        double* GetPressure();
        double GetPressure(int iNode);
        double GetTimeStep(int iNode);
        double* GetTemperature();
        double GetTemperature(int iNode);
        double* GetConservative(int iEqu);
        double GetCons(int iEqu, int iNode);
        double* GetResidual(int iEqu);
        double GetResidual(int iEqu, int iNode);
        double* GetLimiter(int iEqu);
        double GetLimiter(int iEqu, int iNode);
        double* GetViscousFlux(int iEqu, int iDim);
        double GetViscousFlux(int iEqu, int iDim, int iNode);
        double* GetViscousFluxGrad(int iEqu, int iDim, int iGradDim);
        double GetViscousFluxGrad(int iEqu, int iDim, int iGradDim, int iNode);
        double* GetPrimitiveGrad(int iEqu, int iGradDim);
        double GetPrimGrad(int iEqu, int iGradDim, int iNode);
        int* GetNonPhysical();
        int GetNonPhysical(int iNode);
    private:
        int m_equ_num = 5;
        /// @brief 基本变量, 1D 变量索引, 2D 节点索引
        double** m_prim;
        /// @brief 守恒变量, 1D 变量索引, 2D 节点索引
        double** m_cons;
        /// @brief 时间步长
        double* m_dt;
        /// @brief 残差, 1D 变量索引, 2D 节点索引
        double** m_residual;
        /// @brief 限制器, 1D 变量索引, 2D 节点索引
        double** m_limiter;
        /// @brief 温度
        double* m_temperture;
        /// @brief 温度梯度, 1D 梯度维度索引, 2D 节点索引
        double** m_temperture_grad;
        /// @brief 基本变量梯度, 1D 变量索引, 2D 梯度维度索引, 3D 节点索引
        double*** m_prim_grad;
        /// @brief 粘性通量, 1D 变量索引, 2D 通量维度索引, 3D 节点索引
        double*** m_viscous_flux;
        /// @brief 粘性通量梯度, 1D 变量索引, 2D 通量维度索引, 3D 梯度维度索引, 4D 节点索引
        double**** m_viscous_flux_grad;
        /// @brief 非物理点标记
        int* m_non_physical;
    };

}