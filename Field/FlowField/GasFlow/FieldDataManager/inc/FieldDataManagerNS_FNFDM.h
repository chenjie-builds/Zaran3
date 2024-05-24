#pragma once
#include"FieldData.h"
#include "GridFNFDM.h"
#include "FieldDataManager.h"
namespace zaran
{
    class DataManagerNS_FNFDM :public DataManager
    {
    public:
        DataManagerNS_FNFDM(GridFN* grid, FieldData* fieldData, int equ_num);
        ~DataManagerNS_FNFDM();
    public:
        void CreateFieldData()override;
        void RegisterFieldData()override;
    public:
        void SetPrimitive(int iEqu, int iNode, double value);
        void SetConservative(int iEqu, int iNode, double value);
        void SetResidual(int iEqu, int iNode, double value);
        void SetLimiter(int iEqu, int iNode, double value);
        void SetViscousFlux(int iEqu, int iDim, int iNode, double value);
        void SetViscousFluxGrad(int iEqu, int iDim, int iGradDim, int iNode, double value);
        void SetPrimitiveGrad(int iEqu, int iGradDim, int iNode, double value);
        void SetTemperature(int iNode, double value);
        void SetTemperatureGrad(int iGradDim, int iNode, double value);
        void SetNonPhysical(int iNode, int value);
        void SetTimeStep(int iNode, double dt);
    public:
        double GetPrimitive(int iEqu, int iNode);
        double * GetPrimitive(int iEqu);
        double* GetDensity();
        double GetDensity(int iNode);
        double* GetVelocity(int iDim);
        double GetVelocity(int iDim, int iNode);
        double* GetPressure();
        double GetPressure(int iNode);
        double GetTimeStep(int iNode);
        double* GetTemperature();
        double GetTemperature(int iNode);
        double* GetConservative(int iEqu);
        double GetConservative(int iEqu, int iNode);
        double* GetResidual(int iEqu);
        double GetResidual(int iEqu, int iNode);
        double* GetLimiter(int iEqu);
        double GetLimiter(int iEqu, int iNode);
        double* GetViscousFlux(int iEqu, int iDim);
        double GetViscousFlux(int iEqu, int iDim, int iNode);
        double* GetViscousFluxGrad(int iEqu, int iDim, int iGradDim);
        double GetViscousFluxGrad(int iEqu, int iDim, int iGradDim, int iNode);
        double* GetPrimitiveGrad(int iEqu, int iGradDim);
        double GetPrimitiveGrad(int iEqu, int iGradDim, int iNode);
        int* GetNonPhysical();
        int GetNonPhysical(int iNode);

    private:
        GridFN* m_grid;
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