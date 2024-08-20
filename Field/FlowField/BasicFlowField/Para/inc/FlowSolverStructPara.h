#pragma once
#include "FlowSolverPara.h"
namespace zaran
{
    enum class MetricsScheme
    {
        DEER,
        CMM,
    };
    enum class MetricsType
    {
        Originnal,
        CMM1,
        CMM2,
        SCMM,
    };
    enum class FluxDifferenceScheme
    {
        SecondOrder,
        SixthOrder,
    };
    enum class InterSchme
    {
        FirstOrder,
        MUSCL,
        WCNS5,
    };


    class FlowSolverStructPara :public FlowSolverPara
    {
    public:
        FlowSolverStructPara();
        virtual ~FlowSolverStructPara();
        void Init() override;
        void InitMetricsScheme();
        void InitMetricsType();
        void InitFluxDifferenceOrder();
        void InitInterSchme();
    public:
        void SetMetricsScheme(MetricsScheme& metrics_scheme);
        void SetMetricsType(MetricsType& metrics_type);
        void SetFluxDifferenceScheme(FluxDifferenceScheme& flux_difference_scheme);
        void SetInterSchme(InterSchme& inter_schme);
    public:
        const MetricsScheme& GetMetricsScheme()const;
        const MetricsType& GetMetricsType()const;
        const FluxDifferenceScheme& GetFluxDifferenceScheme()const;
        const InterSchme& GetInterSchme()const;
    private:
        MetricsScheme m_metrics_scheme;
        MetricsType m_metrics_type;
        FluxDifferenceScheme m_flux_difference_order;
        InterSchme m_inter_schme;

    };
}