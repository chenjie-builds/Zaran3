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
        InitMetricsScheme();
        InitMetricsType();
        InitFluxDifferenceOrder();
        InitInterSchme();
    }

    void FlowSolverStructPara::InitMetricsScheme()
    {
        string metrics_scheme = GlobalData::GetString("metrics_scheme");
        if (metrics_scheme == "DEER")
        {
            m_metrics_scheme = MetricsScheme::DEER;
        }
        else if (metrics_scheme == "CMM")
        {
            m_metrics_scheme = MetricsScheme::CMM;
        }
        else
        {
            Log::warn("MetricsScheme is not defined, use DEER as default");
            m_metrics_scheme = MetricsScheme::DEER;
        }
    }

    void FlowSolverStructPara::InitMetricsType()
    {
        string metrics_type = GlobalData::GetString("metrics_type");
        if (metrics_type == "Original")
        {
            m_metrics_type = MetricsType::Originnal;
        }
        else if (metrics_type == "CMM1")
        {
            m_metrics_type = MetricsType::CMM1;
        }
        else if (metrics_type == "CMM2")
        {
            m_metrics_type = MetricsType::CMM2;
        }
        else if (metrics_type == "SCMM")
        {
            m_metrics_type = MetricsType::SCMM;
        }
        else
        {
            Log::warn("MetricsType is not defined, use Originnal as default");
            m_metrics_type = MetricsType::Originnal;
        }
    }

    void FlowSolverStructPara::InitFluxDifferenceOrder()
    {
        int flux_scheme_order = GlobalData::GetInt("flux_scheme_order");
        if (flux_scheme_order == 2)
        {
            m_flux_difference_order = FluxDifferenceScheme::SecondOrder;
        }
        else if (flux_scheme_order == 6)
        {
            m_flux_difference_order = FluxDifferenceScheme::SixthOrder;
        }
        else
        {
            Log::warn("FluxDifferenceScheme is not defined, use SecondOrder as default");
            m_flux_difference_order = FluxDifferenceScheme::SecondOrder;
        }
    }

    void FlowSolverStructPara::InitInterSchme()
    {
        string inter_schme = GlobalData::GetString("inter_scheme");
        if (inter_schme == "FirstOrder")
        {
            m_inter_schme = InterSchme::FirstOrder;
        }
        else if (inter_schme == "MUSCL")
        {
            m_inter_schme = InterSchme::MUSCL;
        }
        else if (inter_schme == "WCNS5")
        {
            m_inter_schme = InterSchme::WCNS5;
        }
        else
        {
            Log::warn("InterSchme is not defined, use FirstOrder as default");
            m_inter_schme = InterSchme::FirstOrder;
        }
    }


    const MetricsScheme& FlowSolverStructPara::GetMetricsScheme() const
    {
        return m_metrics_scheme;
    }

    const MetricsType& FlowSolverStructPara::GetMetricsType() const
    {
        return m_metrics_type;
    }

    const FluxDifferenceScheme& FlowSolverStructPara::GetFluxDifferenceScheme() const
    {
        return m_flux_difference_order;
    }

    const InterSchme& FlowSolverStructPara::GetInterSchme() const
    {
        return m_inter_schme;
    }

}