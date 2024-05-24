#include "FieldDataManager.h"
#include "FieldDataManagerNS_FNFDM.h"
namespace zaran
{
    DataManagerNS_FNFDM::DataManagerNS_FNFDM(GridFN* grid, FieldData* fieldData, int equ_num) :DataManager( fieldData, equ_num), m_grid(grid)
    {
        CreateFieldData();
        RegisterFieldData();
    }
    DataManagerNS_FNFDM::~DataManagerNS_FNFDM()
    {
        delete[] m_prim;
        delete[] m_cons;
        delete[] m_residual;
        delete[] m_limiter;
        for (int iEqu = 0; iEqu < m_equ_num; iEqu++)
        {
            delete[] m_viscous_flux[iEqu];
            delete[] m_prim_grad[iEqu];
            for (int iDim = 0; iDim < 3; iDim++)
            {
                delete[] m_viscous_flux_grad[iEqu][iDim];
            }
            delete[] m_viscous_flux_grad[iEqu];
        }
        delete[] m_viscous_flux;
        delete[] m_prim_grad;
        delete[] m_viscous_flux_grad;
        delete[] m_dt;
        delete[] m_temperture;
        delete[] m_temperture_grad;
    }

    void DataManagerNS_FNFDM::SetPrimitive(int iEqu, int iNode, double value)
    {
        m_prim[iEqu][iNode] = value;
    }

    void DataManagerNS_FNFDM::SetConservative(int iEqu, int iNode, double value)
    {
        m_cons[iEqu][iNode] = value;
    }

    void DataManagerNS_FNFDM::SetResidual(int iEqu, int iNode, double value)
    {
        m_residual[iEqu][iNode] = value;
    }

    void DataManagerNS_FNFDM::SetLimiter(int iEqu, int iNode, double value)
    {
        m_limiter[iEqu][iNode] = value;
    }

    void DataManagerNS_FNFDM::SetViscousFlux(int iEqu, int iDim, int iNode, double value)
    {
        m_viscous_flux[iEqu][iDim][iNode] = value;
    }

    void DataManagerNS_FNFDM::SetViscousFluxGrad(int iEqu, int iDim, int iGradDim, int iNode, double value)
    {
        m_viscous_flux_grad[iEqu][iDim][iGradDim][iNode] = value;
    }

    void DataManagerNS_FNFDM::SetPrimitiveGrad(int iEqu, int iGradDim, int iNode, double value)
    {
        m_prim_grad[iEqu][iGradDim][iNode] = value;
    }

    void DataManagerNS_FNFDM::SetTemperature(int iNode, double value)
    {
        m_temperture[iNode] = value;
    }

    void DataManagerNS_FNFDM::SetTemperatureGrad(int iGradDim, int iNode, double value)
    {
        m_temperture_grad[iGradDim][iNode] = value;
    }

    void DataManagerNS_FNFDM::SetNonPhysical(int iNode, int value)
    {
        m_non_physical[iNode] = value;
    }

    void DataManagerNS_FNFDM::SetTimeStep(int iNode, double dt)
    {
        m_dt[iNode] = dt;
    }

    double DataManagerNS_FNFDM::GetPrimitive(int iEqu, int iNode)
    {
        return m_prim[iEqu][iNode];
    }

    double* DataManagerNS_FNFDM::GetPrimitive(int iEqu)
    {
        return  m_prim[iEqu];
    }

    double DataManagerNS_FNFDM::GetDensity(int iNode)
    {
        return m_prim[0][iNode];
    }

     double* DataManagerNS_FNFDM::GetDensity()
    {
        return m_prim[0];
    }

     double* DataManagerNS_FNFDM::GetVelocity(int iDim)
    {
        return m_prim[iDim + 1];
    }

    double DataManagerNS_FNFDM::GetVelocity(int iDim, int iNode)
    {
        return m_prim[iDim + 1][iNode];
    }

     double* DataManagerNS_FNFDM::GetPressure()
    {
        return m_prim[m_equ_num - 1];
    }

    double DataManagerNS_FNFDM::GetPressure(int iNode)
    {
        return m_prim[m_equ_num - 1][iNode];
    }

    double DataManagerNS_FNFDM::GetTimeStep(int iNode)
    {
        return m_dt[iNode];
    }

     double* DataManagerNS_FNFDM::GetTemperature()
    {
        return m_temperture;
    }

    double DataManagerNS_FNFDM::GetTemperature(int iNode)
    {
        return m_temperture[iNode];
    }

     double* DataManagerNS_FNFDM::GetConservative(int iEqu)
    {
        return  m_cons[iEqu];
    }

    double DataManagerNS_FNFDM::GetConservative(int iEqu, int iNode)
    {
        return m_cons[iEqu][iNode];
    }

     double* DataManagerNS_FNFDM::GetResidual(int iEqu)
    {
        return m_residual[iEqu];
    }

    double DataManagerNS_FNFDM::GetResidual(int iEqu, int iNode)
    {
        return m_residual[iEqu][iNode];
    }

     double* DataManagerNS_FNFDM::GetLimiter(int iEqu)
    {
        return m_limiter[iEqu];
    }

    double DataManagerNS_FNFDM::GetLimiter(int iEqu, int iNode)
    {
        return m_limiter[iEqu][iNode];
    }

     double* DataManagerNS_FNFDM::GetViscousFlux(int iEqu, int iDim)
    {
        return m_viscous_flux[iEqu][iDim];
    }

    double DataManagerNS_FNFDM::GetViscousFlux(int iEqu, int iDim, int iNode)
    {
        return m_viscous_flux[iEqu][iDim][iNode];
    }

     double* DataManagerNS_FNFDM::GetViscousFluxGrad(int iEqu, int iDim, int iGradDim)
    {
        return  m_viscous_flux_grad[iEqu][iDim][iGradDim];
    }

    double DataManagerNS_FNFDM::GetViscousFluxGrad(int iEqu, int iDim, int iGradDim, int iNode)
    {
        return m_viscous_flux_grad[iEqu][iDim][iGradDim][iNode];
    }

     double* DataManagerNS_FNFDM::GetPrimitiveGrad(int iEqu, int iGradDim)
    {
        return m_prim_grad[iEqu][iGradDim];
    }

    double DataManagerNS_FNFDM::GetPrimitiveGrad(int iEqu, int iGradDim, int iNode)
    {
        return m_prim_grad[iEqu][iGradDim][iNode];
    }

     int* DataManagerNS_FNFDM::GetNonPhysical()
    {
        return m_non_physical;
    }

    int DataManagerNS_FNFDM::GetNonPhysical(int iNode)
    {
        return m_non_physical[iNode];
    }



    void DataManagerNS_FNFDM::CreateFieldData()
    {
        int node_num = m_grid->GetTotalNodeNum();
        FieldDataType type = FieldDataType::real;
        m_data->AddData("primitive", type, node_num * m_equ_num);
        m_data->AddData("conservative", type, node_num * m_equ_num);
        m_data->AddData("dt", type, node_num);
        m_data->AddData("residual", type, node_num * m_equ_num);
        m_data->AddData("limiter", type, node_num * m_equ_num);
        m_data->AddData("temperature", type, node_num);
        m_data->AddData("temperture_grad", type, node_num * 3);
        m_data->AddData("prim_grad", type, node_num * m_equ_num * 3);
        m_data->AddData("viscous_flux", type, node_num * m_equ_num * 3);
        m_data->AddData("viscous_flux_grad", type, node_num * m_equ_num * 3 * 3);
        m_data->AddData("non_physical", FieldDataType::integer, node_num);
    }
    void DataManagerNS_FNFDM::RegisterFieldData()
    {
        m_prim = new double* [m_equ_num];
        m_cons = new double* [m_equ_num];
        m_residual = new double* [m_equ_num];
        m_limiter = new double* [m_equ_num];
        m_viscous_flux = new double** [m_equ_num];
        m_prim_grad = new double** [m_equ_num];
        m_viscous_flux_grad = new double*** [m_equ_num];
        for (int iEqu = 0; iEqu < m_equ_num; iEqu++)
        {
            m_viscous_flux[iEqu] = new double* [3];
            m_prim_grad[iEqu] = new double* [3];
            m_viscous_flux_grad[iEqu] = new double** [3];
            for (int iDim = 0; iDim < 3; iDim++)
            {
                m_viscous_flux_grad[iEqu][iDim] = new double* [3];
            }
        }
        m_data->GetData("primitive", m_prim[0]);
        m_data->GetData("conservative", m_cons[0]);
        m_data->GetData("residual", m_residual[0]);
        m_data->GetData("limiter", m_limiter[0]);
        for (int iEqu = 1; iEqu < m_equ_num; iEqu++)
        {
            m_prim[iEqu] = m_prim[iEqu - 1] + m_grid->GetTotalNodeNum();
            m_cons[iEqu] = m_cons[iEqu - 1] + m_grid->GetTotalNodeNum();
            m_residual[iEqu] = m_residual[iEqu - 1] + m_grid->GetTotalNodeNum();
            m_limiter[iEqu] = m_limiter[iEqu - 1] + m_grid->GetTotalNodeNum();
        }

        m_data->GetData("viscous_flux", m_viscous_flux[0][0]);
        double* temp_data_ptr;
        for (int iEqu = 0;iEqu < m_equ_num;iEqu++)
        {
            for (int iDim = 0;iDim < 3;iDim++)
            {
                if (iEqu == 0 && iDim == 0)
                    temp_data_ptr = m_viscous_flux[0][0];
                else
                    temp_data_ptr += m_grid->GetTotalNodeNum();
                m_viscous_flux[iEqu][iDim] = temp_data_ptr;
            }
        }
        m_data->GetData("prim_grad", m_prim_grad[0][0]);
        for (int iEqu = 0;iEqu < m_equ_num;iEqu++)
        {
            for (int iDim = 0;iDim < 3;iDim++)
            {
                if (iEqu == 0 && iDim == 0)
                    temp_data_ptr = m_prim_grad[0][0];
                else
                    temp_data_ptr += m_grid->GetTotalNodeNum();
                m_prim_grad[iEqu][iDim] = temp_data_ptr;
            }
        }
        m_data->GetData("viscous_flux_grad", m_viscous_flux_grad[0][0][0]);
        for (int iEqu = 0;iEqu < m_equ_num;iEqu++)
        {
            for (int iDim = 0;iDim < 3;iDim++)
            {
                for (int jDim = 0;jDim < 3;jDim++)
                {
                    if (iEqu == 0 && iDim == 0 && jDim == 0)
                        temp_data_ptr = m_viscous_flux_grad[0][0][0];
                    else
                        temp_data_ptr += m_grid->GetTotalNodeNum();
                    m_viscous_flux_grad[iEqu][iDim][jDim] = temp_data_ptr;
                }
            }
        }
        m_data->GetData("dt", m_dt);
        m_data->GetData("temperture", m_temperture);
        m_data->GetData("non_physical", m_non_physical);
        m_temperture_grad = new double* [3];
        m_data->GetData("temperture_grad", m_temperture_grad[0]);
        for (int iDim = 1;iDim < 3;iDim++)
        {
            m_temperture_grad[iDim] = m_temperture_grad[iDim - 1] + m_grid->GetTotalNodeNum();
        }


    }
}