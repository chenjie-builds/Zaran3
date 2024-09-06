#pragma once
#include "FlowSolverPara.h"
namespace zaran
{
    //半点度量计算方法
    enum class MidMetricsScheme
    {
        DEER,
        CMM,
    };
    //度量计算方法
    //S0:原始格式
    //S1:守恒形式1
    //S2:守恒形式2
    //S3:SCMM格式
    enum class MetricsScheme
    {
        S0,
        S1,
        S2,
        S3,
    };
    //通量差分格式
    //SecondOrder:二阶格式
    //SixthOrder:六阶格式
    enum class DifferenceScheme
    {
        SecondOrder,
        SixthOrder,
    };
    //插值格式
    //FirstOrder:一阶插值
    //MUSCL: MUSCL插值
    //WCNS5:WCNS5插值
    //Grad:梯度插值
    enum class InterpolationScheme
    {
        FirstOrder,
        MUSCL,
        WCNS5,
        Grad,
    };
    //Jacobi计算格式
    //V1:原始格式
    //V2:守恒形式1
    //V3:守恒形式2
    enum class JacobianScheme
    {
        V1,
        V2,
        V3,
    };

    class FlowSolverStructPara :public FlowSolverPara
    {
    public:
        FlowSolverStructPara();
        virtual ~FlowSolverStructPara();
        void Init() override;
        void InitMidMetricsScheme();
        void InitMetricsType();
        void InitFluxDifferenceOrder();
        void InitInterSchme();
        void InitJacobianScheme();
    public:
        const MidMetricsScheme& GetMidMetricsScheme()const;
        const MetricsScheme& GetMetricsScheme()const;
        const DifferenceScheme& GetDifferenceScheme()const;
        const InterpolationScheme& GetInterSchme()const;
        const JacobianScheme& GetJacobianScheme()const;
    private:
        MidMetricsScheme m_mid_metrics_scheme;
        MetricsScheme m_metrics_scheme;
        DifferenceScheme m_difference_scheme;
        InterpolationScheme m_inter_schme;
        JacobianScheme m_jacobian_scheme;
    };
}