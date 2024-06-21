#pragma once
#include "NSSolver.h"
#include "GridStruct.h"
#include"StructIdxProxy.h"
#include "DataManagerNS.h"
#include "NodeMetricsFN.h"
namespace zaran
{
    class NSSolverStruct :public NSSolver
    {
    public:
        NSSolverStruct(int index, string name, FlowSolverPara* para, GridStruct* grid, FieldData* fieldData, DataManagerNS* data_manager);
        ~NSSolverStruct();
    protected:
        void InitFieldFarFlow()override;
        void InitFieldFarFieldNoVelocity()override;
        void InitFieldBackup()override;
        void CalcMetric()override;
    protected:
        void Preprocess()override;
        void CalcTimeStepLocal() override;
        void ReduceTimeStep(double& dt)override;
        void RungeKutta()override;
        void Prim2Cons()override;
        void Cons2Prim()override;
        void ZeroResidual() override;
        /// @brief mid node value reconstruct
        /// @param idx_left left node index
        /// @param idx_right right node index 
        /// @param value_rec_left reconstructed value of left side at mid point
        /// @param value_rec_right reconstructed value of right side at mid point
        void MidNodeRec2nd(int idx_left, int idx_right, double* rec_value_left, double* rec_value_right);
        void MidNodeRec1st(int idx_left, int idx_right, double* rec_value_left, double* rec_value_right);

        // 计算流动通量
        void ConvectiveResidual()override;
        //计算粘性通量
        void ViscousResidual() override;
        void CalcViscousFlux() override;
        void CalcViscousFluxGrad()override;
        //计算源项
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
        DataManagerNS* m_data_manager;
        NodeMetric* m_node_metric;
        StructIdxProxy* m_idx_proxy;

    };
}