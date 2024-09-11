#pragma once
#include "DataManagerNSStruct.h"
#include "FlowSolverStructPara.h"
#include "GridStruct.h"
#include "NSSolver.h"
#include "NodeMetricsFN.h"
#include "StructIdxProxy.h"

namespace zaran
{
  class NSSolverStruct : public NSSolver
  {
  public:
    NSSolverStruct(int index, string name, FlowSolverPara *para, GridStruct *grid, DataManagerNSStruct *data_manager);
    ~NSSolverStruct();
    DataManagerNSStruct *GetDataManager() override;
    Metrics *GetNodeMetrics();
    StructIdxProxy *GetIdxProxy();
    FlowSolverStructPara *GetPara() override;

  protected:
    void InitFieldFarFlow() override;
    void InitFieldFarFieldNoVelocity() override;
    void InitFieldBackup() override;
    void InitFieldVortex();

  protected:
    void CalcCoordTransCoef() override;
    void CalcMetrics();
    /// @brief 使用原始形式计算度量系数
    void CalcMetricsS0();
    void CalcMetricsS0_2nd();
    void CalcMetricsS0_6th();
    /// @brief 计算结构网格节点度量，使用S1
    virtual void CalcMetricsS1();
    void CalcMetricsS1_2nd();
    void CalcMetricsS1_6th();
    /// @brief 计算结构网格节点度量，使用S2
    virtual void CalcMetricsS2();
    void CalcMetricsS2_2nd();
    void CalcMetricsS2_6th();
    /// @brief 计算结构网格节点度量，使用SCMM
    virtual void CalcMetricsS3();
    void CalcMetricsS3_2nd();
    void CalcMetricsS3_6th();
    /// @brief 计算结构网格节点度量，使用FSCMM
    virtual void CalcMetricsS4();
    void CalcJacobian();
    void CalcJacobianV1();
    void CalcJacobianV2();
    void CalcJacobianV3();
    void CalcJacobianV3_2nd();
    void CalcJacobianV3_6th();

  protected:
    void Preprocess() override;
    void CalcPrimGradWLS();
    void CalcLimiterVK();
    void CalcTimeStepLocal() override;
    void ReduceTimeStep(double &dt) override;
    void RungeKutta() override;
    void Prim2Cons() override;
    void Cons2Prim() override;
    void ZeroResidual() override;
    /// @brief 半点原始变量插值：一阶
    void CalcMidNode1st();
    /// @brief 半点原始变量插值：梯度二阶
    void CalcMidNodeGrad();
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
    void MidNodeGrad(int idx_left, int idx_right, double *coord_vec, double *value_left, double *value_right);
    void MidNode1st(int idx_left, int idx_right, double *value_left, double *value_right);
    /// @brief mid node (i+1/2) value interpolate use MUSCL
    /// @param value values at node, size = 5,(i-2,i-1,i,i+1,i+2)
    /// @param value_left interpolated value of left side at mid point i+1/2
    /// @param value_right interpolated value of right side at mid point i+1/2
    void MidNodeMUSCL(const double *value, double &value_left, double &value_right);
    /// @brief mid node value interpolate use WCNS5
    /// @param value values at node, size = 5,(i-2,i-1,i,i+1,i+2)
    /// @param value_left interpolated value of left side at mid point i+1/2
    /// @param value_right interpolated value of right side at mid point i+1/2
    void MidNodeWCNS5(const double *value, double &value_left, double &value_right);
    virtual void FluxDifference2nd() {};
    virtual void FluxDifference6th() {};

    // 计算流动通量
    void CalcInviscidResidual() override;
    virtual void CalcInviscidResidual1st();
    virtual void CalcInviscidResidualGrad();
    virtual void CalcInviscidResidualMUSCL();
    virtual void CalcInviscidResidualWCNS5();
    // 计算粘性通量
    void CalcViscousResidual() override;
    void CalcViscousFlux() override;
    void CalcViscousFluxGrad() override;
    void SourceResidual() override;

    void BoundaryCondition() override;
    void InletBC(BoundStruct &bound);
    void OutletBC(BoundStruct &bound);
    void WallBC(BoundStruct &bound);
    void RiemannBC(BoundStruct &bound);
    void SymmetryBC(BoundStruct &bound);
    void VortexBC(BoundStruct &bound);

  protected:
    void CheckPrimtive() override;
    void CheckResidual() override;
    void FixPrimtive() override;
    void BackupField(std::string &back_folder) override;

  protected:
    GridStruct *GetGrid() override;
    Metrics *GetMidMetricsI();
    Metrics *GetMidMetricsJ();
    Metrics *GetMidMetricsK();

  protected:
    //@brief 节点差分：二阶中心差分
    //@param mid_data 节点数据模板(i-1/2,i+1/2)
    //@return 差分值
    double NodeDifferece2nd(double *mid_data);
    //@brief 节点差分：四阶中心差分
    //@param mid_data 节点数据模板(i-3/2,i-1/2,i+1/2,i+3/2)
    //@return 差分值(i)
    double NodeDifferece4th(double *mid_data);
    //@brief 节点差分: 四阶左差分，用于右侧节点
    //@param mid_data 节点数据模板(i+1/2,i-1/2,i-3/2,i-5/2,i-7/2)
    //@return 差分值（i）
    double NodeDifferece4thLeft(double *mid_data);
    //@brief 节点差分: 四阶右差分，用于左侧节点
    //@param mid_data 节点数据模板(i-1/2,i+1/2,i+3/2,i+5/2,i+7/2)
    //@return 差分值（i）
    double NodeDifferece4thRight(double *mid_data);
    //@brief 节点差分：六阶中心差分
    //@param mid_data 节点数据模板(i-5/2,i-3/2,i-1/2,i+1/2,i+3/2,i+5/2)
    //@return 差分值(i)
    double NodeDifferece6th(double *mid_data);

    //@brief 半点插值：二阶插值
    //@param node_data 半点数据模板(i,i+1)
    //@return 插值值(i+1/2)
    double MidNodeInter2nd(double *node_data);
    //@brief 半点插值：二阶左插值，用于右侧节点
    //@param node_data 半点数据模板(i-1,i)
    //@return 插值值(i-1/2)
    double MidNodeInter2ndLeft(double *node_data);
    //@brief 半点插值：二阶右插值，用于左侧节点
    //@param node_data 半点数据模板(i,i+1)
    //@return 插值值(i-1/2)
    double MidNodeInter2ndRight(double *node_data);

    //@brief 半点插值：四阶插值
    //@param node_data 半点数据模板(i-1,i,i+1,i+2)
    //@return 插值值(i)
    double MidNodeInter4th(double *node_data);
    //@brief 半点插值：四阶左插值，用于右侧节点
    //@param node_data 半点数据模板(i-3,i-2,i-1,i)
    //@return 插值值(i+1/2)
    double MidNodeInter4thLeft1(double *node_data);
    //@brief 半点插值：四阶左插值，用于右侧节点
    //@param node_data 半点数据模板(i-3,i-2,i-1,i)
    //@return 插值值(i-1/2)
    double MidNodeInter4thLeft2(double *node_data);
    //@brief 半点插值：四阶右插值，用于左侧节点
    //@param node_data 半点数据模板(i,i+1,i+2,i+3)
    //@return 插值值(i-1/2)
    double MidNodeInter4thRight1(double *node_data);
    //@brief 半点插值：四阶右插值，用于左侧节点
    //@param node_data 半点数据模板(i,i+1,i+2,i+3)
    //@return 插值值(i+1/2)
    double MidNodeInter4thRight2(double *node_data);
    //@brief 半点插值：六阶插值
    //@param node_data 整点数据模板(i-2,i-1,i,i+1,i+2,i+3)
    //@return 插值值(i/2)
    double MidNodeInter6th(double *node_data);

  private:
    /// @brief 节点度量
    Metrics *m_node_metrics;
    /// @brief 索引代理，用于将结构节点索引转换为场数据索引
    StructIdxProxy *m_idx_proxy;
    Metrics *m_metrics_half_i; // 半点度量系数(i+1/2)
    Metrics *m_metrics_half_j; // 半点度量系数(j+1/2)
    Metrics *m_metrics_half_k; // 半点度量系数(k+1/2)
  };
} // namespace zaran