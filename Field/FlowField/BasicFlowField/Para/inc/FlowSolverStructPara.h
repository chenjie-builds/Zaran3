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
    enum class MetricType
    {
        S0,
        S1,
        S2,
        S3,
    };
    //通量差分格式
    //SecondOrder:二阶格式
    //SixthOrder:六阶格式
    enum class FluxDifferenceScheme
    {
        SecondOrder,
        SixthOrder,
    };
    // 度量系数差分格式
    enum class MetricDifferenceScheme
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
        MUSCL_SV,
        WCNS5,
        Grad,
    };
    //Jacobi计算格式
    //V1:原始格式
    //V2:守恒形式1
    //V3:守恒形式2
    enum class JacobianType
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
        void InitFluxDifferenceScheme();
        void InitMetricDifferenceScheme();
        void InitInterSchme();
        void InitJacobianType();
    public:
        const MidMetricsScheme& GetMidMetricsScheme()const;
        const MetricType& GetMetricsType()const;
        const FluxDifferenceScheme& GetFluxDifferenceScheme()const;
        const MetricDifferenceScheme& GetMetricDifferenceScheme()const;
        const InterpolationScheme& GetInterSchme()const;
        const JacobianType& GetJacobianType()const;
    private:
        MidMetricsScheme m_mid_metrics_scheme;
        MetricType m_metric_type;
        FluxDifferenceScheme m_flux_difference_scheme;
        MetricDifferenceScheme m_metric_difference_scheme;
        InterpolationScheme m_inter_schme;
        JacobianType m_jacobian_type;
    };
}