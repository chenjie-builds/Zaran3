#include "FieldDataManager.h"
#include "DataManagerNS.h"
namespace zaran
{
    DataManagerNS::DataManagerNS(FieldData* fieldData, int data_num) :DataManager(fieldData, data_num)
    {
        CreateData();
        RegisterData();
    }
    DataManagerNS::~DataManagerNS()
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

    void DataManagerNS::SetPrim(int iEqu, int iNode, double value)
    {
        m_prim[iEqu][iNode] = value;
    }

    void DataManagerNS::SetCons(int iEqu, int iNode, double value)
    {
        m_cons[iEqu][iNode] = value;
    }

    void DataManagerNS::SetResidual(int iEqu, int iNode, double value)
    {
        m_residual[iEqu][iNode] = value;
    }

    void DataManagerNS::SetLimiter(int iEqu, int iNode, double value)
    {
        m_limiter[iEqu][iNode] = value;
    }

    void DataManagerNS::SetViscousFlux(int iEqu, int iDim, int iNode, double value)
    {
        m_viscous_flux[iEqu][iDim][iNode] = value;
    }

    void DataManagerNS::SetViscousFluxGrad(int iEqu, int iDim, int iGradDim, int iNode, double value)
    {
        m_viscous_flux_grad[iEqu][iDim][iGradDim][iNode] = value;
    }

    void DataManagerNS::SetPrimitiveGrad(int iEqu, int iGradDim, int iNode, double value)
    {
        m_prim_grad[iEqu][iGradDim][iNode] = value;
    }

    void DataManagerNS::SetTemperature(int iNode, double value)
    {
        m_temperture[iNode] = value;
    }

    void DataManagerNS::SetTemperatureGrad(int iGradDim, int iNode, double value)
    {
        m_temperture_grad[iGradDim][iNode] = value;
    }

    void DataManagerNS::SetNonPhysical(int iNode, int value)
    {
        m_non_physical[iNode] = value;
    }

    void DataManagerNS::SetTimeStep(int iNode, double dt)
    {
        m_dt[iNode] = dt;
    }

    double DataManagerNS::GetPrim(int iEqu, int iNode)
    {
        return m_prim[iEqu][iNode];
    }

    double* DataManagerNS::GetPrim(int iEqu)
    {
        return  m_prim[iEqu];
    }

    double DataManagerNS::GetDensity(int iNode)
    {
        return m_prim[0][iNode];
    }

    double* DataManagerNS::GetDensity()
    {
        return m_prim[0];
    }

    double* DataManagerNS::GetVelocity(int iDim)
    {
        return m_prim[iDim + 1];
    }

    double DataManagerNS::GetVelocity(int iDim, int iNode)
    {
        return m_prim[iDim + 1][iNode];
    }

    double* DataManagerNS::GetPressure()
    {
        return m_prim[m_equ_num - 1];
    }

    double DataManagerNS::GetPressure(int iNode)
    {
        return m_prim[m_equ_num - 1][iNode];
    }

    double DataManagerNS::GetTimeStep(int iNode)
    {
        return m_dt[iNode];
    }

    double* DataManagerNS::GetTemperature()
    {
        return m_temperture;
    }

    double DataManagerNS::GetTemperature(int iNode)
    {
        return m_temperture[iNode];
    }

    double* DataManagerNS::GetConservative(int iEqu)
    {
        return  m_cons[iEqu];
    }

    double DataManagerNS::GetCons(int iEqu, int iNode)
    {
        return m_cons[iEqu][iNode];
    }

    double* DataManagerNS::GetResidual(int iEqu)
    {
        return m_residual[iEqu];
    }

    double DataManagerNS::GetResidual(int iEqu, int iNode)
    {
        return m_residual[iEqu][iNode];
    }

    double* DataManagerNS::GetLimiter(int iEqu)
    {
        return m_limiter[iEqu];
    }

    double DataManagerNS::GetLimiter(int iEqu, int iNode)
    {
        return m_limiter[iEqu][iNode];
    }

    double* DataManagerNS::GetViscousFlux(int iEqu, int iDim)
    {
        return m_viscous_flux[iEqu][iDim];
    }

    double DataManagerNS::GetViscousFlux(int iEqu, int iDim, int iNode)
    {
        return m_viscous_flux[iEqu][iDim][iNode];
    }

    double* DataManagerNS::GetViscousFluxGrad(int iEqu, int iDim, int iGradDim)
    {
        return  m_viscous_flux_grad[iEqu][iDim][iGradDim];
    }

    double DataManagerNS::GetViscousFluxGrad(int iEqu, int iDim, int iGradDim, int iNode)
    {
        return m_viscous_flux_grad[iEqu][iDim][iGradDim][iNode];
    }

    double* DataManagerNS::GetPrimitiveGrad(int iEqu, int iGradDim)
    {
        return m_prim_grad[iEqu][iGradDim];
    }

    double DataManagerNS::GetPrimGrad(int iEqu, int iGradDim, int iNode)
    {
        return m_prim_grad[iEqu][iGradDim][iNode];
    }

    int* DataManagerNS::GetNonPhysical()
    {
        return m_non_physical;
    }

    int DataManagerNS::GetNonPhysical(int iNode)
    {
        return m_non_physical[iNode];
    }



    void DataManagerNS::CreateData()
    {
        FieldDataType type = FieldDataType::real;
        m_data->AddData("density", type, m_data_num);
        m_data->AddData("x_velocity", type, m_data_num);
        m_data->AddData("y_velocity", type, m_data_num);
        m_data->AddData("z_velocity", type, m_data_num);
        m_data->AddData("pressure", type, m_data_num);
        for (int iEqu = 0;iEqu < 5;iEqu++)
        {
            m_data->AddData("primitive_" + std::to_string(iEqu), type, m_data_num);
            m_data->AddData("conservative_" + std::to_string(iEqu), type, m_data_num);
            m_data->AddData("residual_" + std::to_string(iEqu), type, m_data_num);
            m_data->AddData("limiter_" + std::to_string(iEqu), type, m_data_num);
            for (int iDim = 0;iDim < 3;iDim++)
            {
                m_data->AddData("viscous_flux_" + std::to_string(iEqu) + "_" + std::to_string(iDim), type, m_data_num);
                m_data->AddData("prim_grad_" + std::to_string(iEqu) + "_" + std::to_string(iDim), type, m_data_num);
                for (int jDim = 0;jDim < 3;jDim++)
                {
                    m_data->AddData("viscous_flux_grad_" + std::to_string(iEqu) + "_" + std::to_string(iDim) + "_" + std::to_string(jDim), type, m_data_num);
                }
            }
        }
        for (int iDim = 0;iDim < 3;iDim++)
        {
            m_data->AddData("temperture_grad_" + std::to_string(iDim), type, m_data_num);
        }
        m_data->AddData("dt", type, m_data_num);
        m_data->AddData("temperature", type, m_data_num);
        m_data->AddData("non_physical", FieldDataType::integer, m_data_num);
    }
    void DataManagerNS::RegisterData()
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
        m_data->GetData("dt", m_dt);
        m_data->GetData("temperture", m_temperture);
        m_data->GetData("non_physical", m_non_physical);
        for (int iEqu = 0;iEqu < 5;iEqu++)
        {
            m_data->GetData("primitive_" + std::to_string(iEqu), m_prim[iEqu]);
            m_data->GetData("conservative_" + std::to_string(iEqu), m_cons[iEqu]);
            m_data->GetData("residual_" + std::to_string(iEqu), m_residual[iEqu]);
            m_data->GetData("limiter_" + std::to_string(iEqu), m_limiter[iEqu]);
            for (int iDim = 0;iDim < 3;iDim++)
            {
                m_data->GetData("viscous_flux_" + std::to_string(iEqu) + "_" + std::to_string(iDim), m_viscous_flux[iEqu][iDim]);
                m_data->GetData("prim_grad_" + std::to_string(iEqu) + "_" + std::to_string(iDim), m_prim_grad[iEqu][iDim]);
                for (int jDim = 0;jDim < 3;jDim++)
                {
                    m_data->GetData("viscous_flux_grad_" + std::to_string(iEqu) + "_" + std::to_string(iDim) + "_" + std::to_string(jDim), m_viscous_flux_grad[iEqu][iDim][jDim]);
                }
            }
        }
        m_temperture_grad = new double* [3];
        for (int iDim = 0;iDim < 3;iDim++)
        {
            m_data->GetData("temperture_grad_" + std::to_string(iDim), m_temperture_grad[iDim]);
        }



    }
}