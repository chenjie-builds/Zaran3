#include "FlowSolverStructPara.h"
#include "GlobalData.h"
#include "Log.h"
namespace zaran
{



    FlowSolverStructPara::FlowSolverStructPara()
    {
    }

    FlowSolverStructPara::~FlowSolverStructPara()
    {
    }

    void FlowSolverStructPara::Init()
    {
        FlowSolverPara::Init();
        InitMidMetricsScheme();
        InitMetricsType();
        InitFluxDifferenceOrder();
        InitInterSchme();
    }

    void FlowSolverStructPara::InitMidMetricsScheme()
    {
        string mid_metric_method = GlobalData::GetString("mid_metric_method");
        if (mid_metric_method == "DEER")
        {
            m_mid_metrics_scheme = MidMetricsScheme::DEER;
        }
        else if (mid_metric_method == "CMM")
        {
            m_mid_metrics_scheme = MidMetricsScheme::CMM;
        }
        else
        {
            Log::warn("MetricsScheme is not defined, use DEER as default");
            m_mid_metrics_scheme = MidMetricsScheme::DEER;
        }
    }

    void FlowSolverStructPara::InitMetricsType()
    {
        string metrics_scheme = GlobalData::GetString("metrics_scheme");
        if (metrics_scheme == "Original")
        {
            m_metrics_scheme = MetricsScheme::S0;
        }
        else if (metrics_scheme == "CMM1")
        {
            m_metrics_scheme = MetricsScheme::S1;
        }
        else if (metrics_scheme == "CMM2")
        {
            m_metrics_scheme = MetricsScheme::S2;
        }
        else if (metrics_scheme == "SCMM")
        {
            m_metrics_scheme = MetricsScheme::S3;
        }
        else
        {
            Log::warn("MetricsType is not defined, use Originnal as default");
            m_metrics_scheme = MetricsScheme::S0;
        }
    }

    void FlowSolverStructPara::InitFluxDifferenceOrder()
    {
        int difference_scheme = GlobalData::GetInt("difference_scheme");
        if (difference_scheme == 2)
        {
            m_difference_scheme = DifferenceScheme::SecondOrder;
        }
        else if (difference_scheme == 6)
        {
            m_difference_scheme = DifferenceScheme::SixthOrder;
        }
        else
        {
            Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
            m_difference_scheme = DifferenceScheme::SecondOrder;
        }
    }

    void FlowSolverStructPara::InitInterSchme()
    {
        string inter_schme = GlobalData::GetString("inter_scheme");
        if (inter_schme == "FirstOrder")
        {
            m_inter_schme = InterpolationScheme::FirstOrder;
        }
        else if (inter_schme == "MUSCL")
        {
            m_inter_schme = InterpolationScheme::MUSCL;
        }
        else if (inter_schme == "WCNS5")
        {
            m_inter_schme = InterpolationScheme::WCNS5;
        }
        else
        {
            Log::warn("InterSchme is not defined, use FirstOrder as default");
            m_inter_schme = InterpolationScheme::FirstOrder;
        }
    }

    void FlowSolverStructPara::InitJacobianScheme()
    {
        string jacobian_scheme = GlobalData::GetString("jacobian_scheme");
        if (jacobian_scheme == "V1")
        {
            m_jacobian_scheme = JacobianScheme::V1;
        }
        else if (jacobian_scheme == "V2")
        {
            m_jacobian_scheme = JacobianScheme::V2;
        }
        else if (jacobian_scheme == "V3")
        {
            m_jacobian_scheme = JacobianScheme::V3;
        }
        else
        {
            Log::warn("JacobianScheme is not defined, use V1 as default");
            m_jacobian_scheme = JacobianScheme::V1;
        }
    }

    const MidMetricsScheme& FlowSolverStructPara::GetMidMetricsScheme() const
    {
        return m_mid_metrics_scheme;
    }

    const MetricsScheme& FlowSolverStructPara::GetMetricsScheme() const
    {
        return m_metrics_scheme;
    }

    const DifferenceScheme& FlowSolverStructPara::GetDifferenceScheme() const
    {
        return m_difference_scheme;
    }

    const InterpolationScheme& FlowSolverStructPara::GetInterSchme() const
    {
        return m_inter_schme;
    }

    const JacobianScheme &FlowSolverStructPara::GetJacobianScheme() const
    {
        return m_jacobian_scheme;
    }
}