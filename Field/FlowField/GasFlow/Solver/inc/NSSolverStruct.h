#pragma once
#include "NSSolver.h"
#include "GridStruct.h"
#include"StructIdxProxy.h"
#include "DataManagerNSStruct.h"
#include "NodeMetricsFN.h"
namespace zaran
{
    class NSSolverStruct :public NSSolver
    {
    public:
        NSSolverStruct(int index, string name, FlowSolverPara* para, GridStruct* grid, DataManagerNSStruct* data_manager);
        ~NSSolverStruct();
        DataManagerNSStruct* GetDataManager()override;
    protected:
        void InitFieldFarFlow()override;
        void InitFieldFarFieldNoVelocity()override;
        void InitFieldBackup()override;
        void CalcMetrics()override;
        /// @brief 计算结构网格节点度量，使用CMM1
        void CalcMetricsCMM1();
        /// @brief 计算结构网格节点度量，使用CMM2
        void CalcMetricsCMM2();
        /// @brief 计算结构网格节点度量，使用SCMM
        void CalcMetricsSCMM();
        /// @brief 计算结构网格节点度量，使用FSCMM
        void CalcMetricsFSCMM();
    protected:
        void Preprocess()override;
        void CalcTimeStepLocal() override;
        void ReduceTimeStep(double& dt)override;
        void RungeKutta()override;
        void Prim2Cons()override;
        void Cons2Prim()override;
        void ZeroResidual() override;
        /// @brief 半点原始变量插值
        void CalcMidNodePrim();
        /// @brief 半点原始变量插值：一阶
        void CalcMidNode1st();
        /// @brief 半点原始变量插值：梯度二阶
        void CalcMidNodePrimGrad();
        /// @brief 半点原始变量插值：MUSCL
        void CalcMidNodePrimMUSCL();
        /// @brief 虚拟节点半点原始变量插值：MUSCL
        void CalcMidGhostNodePrimMUSCL();
        /// @brief 半点原始变量插值：WCNS5
        void CalcMidNodePrimWCNS5();
        /// @brief 虚拟节点半点原始变量插值：WCNS5
        void CalcMidGhostNodePrimWCNS5();
        /// @brief mid node value interpolate
        /// @param idx_left left node index
        /// @param idx_right right node index 
        /// @param coord_vec coordinate vector from left to right
        /// @param value_left interpolated value of left side at mid point
        /// @param value_right interpolated value of right side at mid point
        void MidNodeGrad(int idx_left, int idx_right, double* coord_vec, double* value_left, double* value_right);
        void MidNode1st(int idx_left, int idx_right, double* value_left, double* value_right);
        /// @brief mid node (i+1/2) value interpolate use MUSCL
        /// @param value values at node, size = 5,(i-2,i-1,i,i+1,i+2)
        /// @param value_left interpolated value of left side at mid point i+1/2
        /// @param value_right interpolated value of right side at mid point i+1/2
        void MidNodeMUSCL(const double* value, double& value_left, double& value_right);
        /// @brief mid node value interpolate use WCNS5
        /// @param value values at node, size = 5,(i-2,i-1,i,i+1,i+2)
        /// @param value_left interpolated value of left side at mid point i+1/2
        /// @param value_right interpolated value of right side at mid point i+1/2
        void MidNodeWCNS5(const double* value, double& value_left, double& value_right);

        void FluxDifferenceMUSCL();
        void FluxDifferenceWCNS5();

        // 计算流动通量
        void ConvectiveResidual()override;
        void ConvectiveResidual1st();
        void ConvectiveResidualGrad();
        void ConvectiveResidualMUSCL();
        void ConvectiveResidualWCNS5();
        //计算粘性通量
        void ViscousResidual() override;
        void CalcViscousFlux() override;
        void CalcViscousFluxGrad()override;
        void SourceResidual() override;

        void BoundaryCondition()override;
        void InletBC(BoundStruct& bound);
        void OutletBC(BoundStruct& bound);
        void WallBC(BoundStruct& bound);
        void RiemannBC(BoundStruct& bound);
        void SymmetryBC(BoundStruct& bound);
    protected:
        void CheckPrimtive()override;
        void CheckResidual()override;
        void FixPrimtive()override;
        void BackupField(std::string& back_folder)override;
    protected:
        virtual GridStruct* GetGrid();
    private:
        /// @brief 节点度量
        Metrics* m_node_metrics;
        /// @brief 索引代理，用于将结构节点索引转换为场数据索引
        StructIdxProxy* m_idx_proxy;
    };
}