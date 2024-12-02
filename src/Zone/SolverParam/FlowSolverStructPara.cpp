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
        InitFluxDifferenceScheme();
        InitMetricDifferenceScheme();
        InitInterScheme();
        InitJacobianType();
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
            Log::warn("MetricType is not defined, use DEER as default");
            m_mid_metrics_scheme = MidMetricsScheme::DEER;
        }
    }

    void FlowSolverStructPara::InitMetricsType()
    {
        string metrics_scheme = GlobalData::GetString("metrics_scheme");
        if (metrics_scheme == "S0")
        {
            m_metric_type = MetricType::S0;
        }
        else if (metrics_scheme == "S1")
        {
            m_metric_type = MetricType::S1;
        }
        else if (metrics_scheme == "S2")
        {
            m_metric_type = MetricType::S2;
        }
        else if (metrics_scheme == "S3")
        {
            m_metric_type = MetricType::S3;
        }
        else
        {
            Log::warn("MetricsType is not defined, use Originnal as default");
            m_metric_type = MetricType::S0;
        }
    }

    void FlowSolverStructPara::InitFluxDifferenceScheme()
    {
        int difference_scheme = GlobalData::GetInt("flux_difference_scheme");
        if (difference_scheme == 2)
        {
            m_flux_difference_scheme = FluxDifferenceScheme::SecondOrder;
        }
		else if (difference_scheme == 4)
		{
			m_flux_difference_scheme = FluxDifferenceScheme::FourthOrder;
		}
        else if (difference_scheme == 6)
        {
            m_flux_difference_scheme = FluxDifferenceScheme::SixthOrder;
        }
        else
        {
            Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
            m_flux_difference_scheme = FluxDifferenceScheme::SecondOrder;
        }
    }

    void FlowSolverStructPara::InitMetricDifferenceScheme()
    {
        int metric_difference_scheme = GlobalData::GetInt("metric_difference_scheme");
        if (metric_difference_scheme == 2)
        {
            m_metric_difference_scheme = MetricDifferenceScheme::SecondOrder;
        }
        else if (metric_difference_scheme == 6)
        {
            m_metric_difference_scheme = MetricDifferenceScheme::SixthOrder;
        }
        else
        {
            Log::warn("MetricDifferenceScheme is not defined, use SecondOrder as default");
            m_metric_difference_scheme = MetricDifferenceScheme::SecondOrder;
        }
    }

    void FlowSolverStructPara::InitInterScheme()
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
        else if (inter_schme == "MUSCL_SV")
        {
            m_inter_schme = InterpolationScheme::MUSCL_SV;
        }
        else if (inter_schme == "WCNS5")
        {
            m_inter_schme = InterpolationScheme::WCNS5;
        }
        else if(inter_schme == "Grad")
        {
            m_inter_schme = InterpolationScheme::Grad;
        }
        else
        {
            Log::warn("Interpolation Scheme:{} is not defined, use First Order as default", inter_schme);
            m_inter_schme = InterpolationScheme::FirstOrder;
        }
    }

    void FlowSolverStructPara::InitJacobianType()
    {
        string jacobian_scheme = GlobalData::GetString("jacobian_scheme");
        if (jacobian_scheme == "V1")
        {
            m_jacobian_type = JacobianType::V1;
        }
        else if (jacobian_scheme == "V2")
        {
            m_jacobian_type = JacobianType::V2;
        }
        else if (jacobian_scheme == "V3")
        {
            m_jacobian_type = JacobianType::V3;
        }
        else
        {
            Log::warn("JacobianType:{} is not defined, use V1 as default", jacobian_scheme);
            m_jacobian_type = JacobianType::V1;
        }
    }

    const MidMetricsScheme& FlowSolverStructPara::GetMidMetricsScheme() const
    {
        return m_mid_metrics_scheme;
    }

    const MetricType& FlowSolverStructPara::GetMetricsType() const
    {
        return m_metric_type;
    }

    const FluxDifferenceScheme& FlowSolverStructPara::GetFluxDifferenceScheme() const
    {
        return m_flux_difference_scheme;
    }

    const MetricDifferenceScheme& FlowSolverStructPara::GetMetricDifferenceScheme() const
    {
        return m_metric_difference_scheme;
    }

    const InterpolationScheme& FlowSolverStructPara::GetInterSchme() const
    {
        return m_inter_schme;
    }

    const JacobianType &FlowSolverStructPara::GetJacobianType() const
    {
        return m_jacobian_type;
    }
}