#pragma once
#include "DataManagerNSStruct.h"
#include "FlowSolverStructPara.h"
#include "GridStruct.h"
#include "NSSolver.h"
#include "Metric.h"
#include "StructIdxProxy.h"

namespace zaran
{
  class NSSolverStruct : public NSSolver
  {
  public:
    NSSolverStruct(int index, string name, FlowSolverParam *para, GridStruct *grid, DataManagerNSStruct *data_manager);
    ~NSSolverStruct();
    DataManagerNSStruct *GetDataManager() override;
    Metric *GetNodeMetrics();
    IdxStruct *GetIdxProxy();
    FlowSolverParamStruct *GetPara() override;

  protected:
    void InitField() override;
    void InitFieldFarfield() override;
    void InitFieldFarFieldZeroVel() override;
    void InitFieldBackup() override;
    void InitFieldVortex();
    void InitFieldExplosion() override;

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
    void CalcPrimGradWLS2D();
    void CalcPrimGradWLS3D();
    void CalcLimiterVK();
    void CalcTimeStepLocal() override;
	void CalcMinTimeStep(double& dt) override;
    void ReduceTimeStep(double &dt) override;
    void RungeKutta() override;
    void Prim2Cons() override;
    void Cons2Prim() override;
    void ZeroResidual() override;
    /// @brief 半点原始变量插值：一阶
    void InterMidNodePrim_1st();
    /// @brief 半点原始变量插值：梯度二阶
    void InterMidNodePrim_Grad();
    /// @brief 半点原始变量插值：MUSCL
    void InterMidNodePrim_MUSCL();
    void InterMidNodePrim_MUSCLSV();
    /// @brief 虚拟节点半点原始变量插值：MUSCL
    void InterGhostMidNode_MUSCL();
    void CalcPrimGhostMidNodeMUSCL_SV();
    /// @brief 半点原始变量插值：WCNS5
    void InterMidNodePrim_WCNS5();
    /// @brief 虚拟节点半点原始变量插值：WCNS5
    void InterGhostMidNodePrim_WCNS5();
    /// @brief mid node primitive value interpolate use gradient
    /// @param idx_left left node index
    /// @param idx_right right node index
    /// @param coord_vec coordinate vector from left to right
    /// @param value_left interpolated value of left side at mid point
    /// @param value_right interpolated value of right side at mid point
    void InterMidNodePrim_Grad(int idx_left, int idx_right, double *coord_vec, double *value_left, double *value_right);
    /// @brief mid node primitive value interpolate use 1st order
    /// @param idx_left 
    /// @param idx_right 
    /// @param value_left 
    /// @param value_right 
    void InterMidNodePrim_1st(int idx_left, int idx_right, double *value_left, double *value_right);
    /// @brief mid node primitive value interpolate use MUSCL
    /// @param value values at node, size = 5,(i-2,i-1,i,i+1,i+2)
    /// @param value_left interpolated value of left side at mid point i+1/2
    /// @param value_right interpolated value of right side at mid point i+1/2
    void InterMidNodePrim_MUSCL(const double *value, double &value_left, double &value_right);
    void InterMidNodePrim_MUSCLSV(const double *value,  double**coord, double &value_left, double &value_right);
    /// @brief mid node primitive value interpolate use WCNS5
    /// @param value values at node, size = 5,(i-2,i-1,i,i+1,i+2)
    /// @param value_left interpolated value of left side at mid point i+1/2
    /// @param value_right interpolated value of right side at mid point i-1/2
    void InterMidNodePrim_WCNS5( double *value, double &value_left, double &value_right);
    void FluxDifference();
    /// @brief flux difference: central 2nd order
    virtual void FluxDifference2nd() {};
	/// @brief flux difference: central 4th order
	virtual void FluxDifference4th() {};
    /// @brief flux difference: central 6th order
    virtual void FluxDifference6th() {};
    /// @brief calculate convection residual
    void CalcConvectionResidual() override;

    virtual void CalcConvectionRes_1st();
    virtual void CalcConvectionRes_Grad();
    virtual void CalcConvectionRes_MUSCL();
    virtual void CalcConvectionRes_MUSCLSV();
    virtual void CalcConvectionRes_WCNS5();

    // 计算粘性通量
    void CalcViscousResidual() override;
    void CalcViscousFlux() override;
    void CalcViscousFluxGrad() override;
    void CalcSourceResidual() override;

    void BoundaryCondition() override;
    void BCInlow(BoundStruct &bound);
    void BCOutflow(BoundStruct &bound);
    void BCWall(BoundStruct &bound);
    void BCFarfield(BoundStruct &bound);
    void BCSymmetry(BoundStruct &bound);
    void BCVortex(BoundStruct &bound);
    //双马赫反射
	void BCDoubleMach(BoundStruct& bound);

  protected:
    void CheckPrimtive() override;
    void CheckResidual() override;
    void FixPrimtive() override;
    void BackupField(std::string &back_folder) override;

  protected:
    GridStruct *GetGrid() override;
    Metric *GetMidMetricsI();
    Metric *GetMidMetricsJ();
    Metric *GetMidMetricsK();

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
    Metric *m_node_metrics;
    /// @brief 索引代理，用于将结构节点索引转换为场数据索引
    IdxStruct *m_idx_proxy;
    Metric *m_metrics_mid_i; // 半点度量系数(i+1/2)
    Metric *m_metrics_mid_j; // 半点度量系数(j+1/2)
    Metric *m_metrics_mid_k; // 半点度量系数(k+1/2)
  };
} // namespace zaran