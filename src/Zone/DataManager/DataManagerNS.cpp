#include "FieldDataManager.h"
#include "DataManagerNS.h"
namespace zaran
{
	DataManagerNS::DataManagerNS(shared_ptr<FieldData> field_data, int data_num) :DataManager(field_data, data_num)
	{
	}
	DataManagerNS::~DataManagerNS()
	{
		delete[] m_prim;
		delete[] m_cons;
		delete[] m_cons_old;
		delete[] m_residual;
		delete[] m_limiter;
		for (int idx_eq = 0; idx_eq < m_equ_num; idx_eq++)
		{
			delete[] m_viscous_flux[idx_eq];
			delete[] m_prim_grad[idx_eq];
			for (int idx_dim = 0; idx_dim < 3; idx_dim++)
			{
				delete[] m_viscous_flux_grad[idx_eq][idx_dim];
			}
			delete[] m_viscous_flux_grad[idx_eq];
		}
		delete[] m_viscous_flux;
		delete[] m_prim_grad;
		delete[] m_viscous_flux_grad;
		delete[] m_dt;
		delete[] m_temperture;
		delete[] m_temperture_grad;
	}

	void DataManagerNS::SetPrim(int idx_eq, int idx_data, double data_value)
	{
		m_prim[idx_eq][idx_data] = data_value;
	}
	void DataManagerNS::SetPrim(int idx_data, const double* value)
	{
		for (int idx_eq = 0; idx_eq < m_equ_num; idx_eq++)
		{
			m_prim[idx_eq][idx_data] = value[idx_eq];
		}
	}
	void DataManagerNS::SetCons(int idx_eq, int idx_data, double data_value)
	{
		m_cons[idx_eq][idx_data] = data_value;
	}
	void DataManagerNS::SetCons(int idx_data, const double* data_value)
	{
		for (int idx_eq = 0; idx_eq < m_equ_num; idx_eq++)
		{
			m_cons[idx_eq][idx_data] = data_value[idx_eq];
		}
	}
	void DataManagerNS::SetConsOld(int idx_eq, int idx_data, double data_value)
	{
		m_cons_old[idx_eq][idx_data] = data_value;
	}

	void DataManagerNS::SetResidual(int idx_eq, int idx_data, double data_value)
	{
		m_residual[idx_eq][idx_data] = data_value;
	}
	void DataManagerNS::SetResidual(int idx_data, const double* value)
	{
		for (int idx_eq = 0; idx_eq < m_equ_num; idx_eq++)
		{
			m_residual[idx_eq][idx_data] = value[idx_eq];
		}
	}
	void DataManagerNS::SetLimiter(int idx_eq, int idx_data, double data_value)
	{
		m_limiter[idx_eq][idx_data] = data_value;
	}

	void DataManagerNS::SetViscousFlux(int idx_eq, int idx_dim, int idx_data, double data_value)
	{
		m_viscous_flux[idx_eq][idx_dim][idx_data] = data_value;
	}

	void DataManagerNS::SetViscousFluxGrad(int idx_eq, int idx_dim, int iGradDim, int idx_data, double data_value)
	{
		m_viscous_flux_grad[idx_eq][idx_dim][iGradDim][idx_data] = data_value;
	}

	void DataManagerNS::SetPrimitiveGrad(int idx_eq, int iGradDim, int idx_data, double data_value)
	{
		m_prim_grad[idx_eq][iGradDim][idx_data] = data_value;
	}

	void DataManagerNS::SetTemperature(int idx_data, double data_value)
	{
		m_temperture[idx_data] = data_value;
	}

	void DataManagerNS::SetTemperatureGrad(int iGradDim, int idx_data, double data_value)
	{
		m_temperture_grad[iGradDim][idx_data] = data_value;
	}

	void DataManagerNS::SetNonPhysical(int idx_data, int data_value)
	{
		m_non_physical[idx_data] = data_value;
	}

	void DataManagerNS::SetTimeStep(int idx_data, double dt)
	{
		m_dt[idx_data] = dt;
	}

	double DataManagerNS::GetPrim(int idx_eq, int idx_data)
	{
		return m_prim[idx_eq][idx_data];
	}

	double* DataManagerNS::GetPrim(int idx_eq)
	{
		return  m_prim[idx_eq];
	}

	
	double DataManagerNS::GetTimeStep(int idx_data)
	{
		return m_dt[idx_data];
	}

	double* DataManagerNS::GetTemperature()
	{
		return m_temperture;
	}

	double DataManagerNS::GetTemperature(int idx_data)
	{
		return m_temperture[idx_data];
	}

	double* DataManagerNS::GetConservative(int idx_eq)
	{
		return  m_cons[idx_eq];
	}
	double* DataManagerNS::GetConsOld(int idx_eq)
	{
		return m_cons_old[idx_eq];
	}

	double DataManagerNS::GetCons(int idx_eq, int idx_data)
	{
		return m_cons[idx_eq][idx_data];
	}
	double DataManagerNS::GetConsOld(int idx_eq, int idx_data)
	{
		return m_cons_old[idx_eq][idx_data];
	}

	double* DataManagerNS::GetResidual(int idx_eq)
	{
		return m_residual[idx_eq];
	}

	double DataManagerNS::GetResidual(int idx_eq, int idx_data)
	{
		return m_residual[idx_eq][idx_data];
	}

	double* DataManagerNS::GetLimiter(int idx_eq)
	{
		return m_limiter[idx_eq];
	}

	double DataManagerNS::GetLimiter(int idx_eq, int idx_data)
	{
		return m_limiter[idx_eq][idx_data];
	}

	double* DataManagerNS::GetViscousFlux(int idx_eq, int idx_dim)
	{
		return m_viscous_flux[idx_eq][idx_dim];
	}

	double DataManagerNS::GetViscousFlux(int idx_eq, int idx_dim, int idx_data)
	{
		return m_viscous_flux[idx_eq][idx_dim][idx_data];
	}

	double* DataManagerNS::GetViscousFluxGrad(int idx_eq, int idx_dim, int iGradDim)
	{
		return  m_viscous_flux_grad[idx_eq][idx_dim][iGradDim];
	}

	double DataManagerNS::GetViscousFluxGrad(int idx_eq, int idx_dim, int iGradDim, int idx_data)
	{
		return m_viscous_flux_grad[idx_eq][idx_dim][iGradDim][idx_data];
	}

	double* DataManagerNS::GetPrimitiveGrad(int idx_eq, int iGradDim)
	{
		return m_prim_grad[idx_eq][iGradDim];
	}

	double DataManagerNS::GetPrimGrad(int idx_eq, int iGradDim, int idx_data)
	{
		return m_prim_grad[idx_eq][iGradDim][idx_data];
	}

	int* DataManagerNS::GetNonPhysical()
	{
		return m_non_physical;
	}

	int DataManagerNS::GetNonPhysical(int idx_data)
	{
		return m_non_physical[idx_data];
	}



	void DataManagerNS::CreateData()
	{
		FieldDataType type = FieldDataType::real;
		// m_data->AddData("density", type, m_data_num);
		// m_data->AddData("x_velocity", type, m_data_num);
		// m_data->AddData("y_velocity", type, m_data_num);
		// m_data->AddData("z_velocity", type, m_data_num);
		// m_data->AddData("pressure", type, m_data_num);
		for (int idx_eq = 0; idx_eq < GetEqNum(); idx_eq++)
		{
			m_data->AddData("primitive_" + std::to_string(idx_eq), type, m_data_num);
			m_data->AddData("conservative_" + std::to_string(idx_eq), type, m_data_num);
			m_data->AddData("conservative_old_" + std::to_string(idx_eq), type, m_data_num);
			m_data->AddData("residual_" + std::to_string(idx_eq), type, m_data_num);
			m_data->AddData("limiter_" + std::to_string(idx_eq), type, m_data_num);
			for (int idx_dim = 0; idx_dim < 3; idx_dim++)
			{
				m_data->AddData("viscous_flux_" + std::to_string(idx_eq) + "_" + std::to_string(idx_dim), type, m_data_num);
				m_data->AddData("prim_grad_" + std::to_string(idx_eq) + "_" + std::to_string(idx_dim), type, m_data_num);
				for (int jDim = 0; jDim < 3; jDim++)
				{
					m_data->AddData("viscous_flux_grad_" + std::to_string(idx_eq) + "_" + std::to_string(idx_dim) + "_" + std::to_string(jDim), type, m_data_num);
				}
			}
		}
		for (int idx_dim = 0; idx_dim < 3; idx_dim++)
		{
			m_data->AddData("temperture_grad_" + std::to_string(idx_dim), type, m_data_num);
		}
		m_data->AddData("dt", type, m_data_num);
		m_data->AddData("temperature", type, m_data_num);
		m_data->AddData("non_physical", FieldDataType::integer, m_data_num);
	}
	void DataManagerNS::RegisterData()
	{
		m_prim = new double* [m_equ_num];
		m_cons = new double* [m_equ_num];
		m_cons_old = new double* [m_equ_num];
		m_residual = new double* [m_equ_num];
		m_limiter = new double* [m_equ_num];
		m_viscous_flux = new double** [m_equ_num];
		m_prim_grad = new double** [m_equ_num];
		m_viscous_flux_grad = new double*** [m_equ_num];
		for (int idx_eq = 0; idx_eq < m_equ_num; idx_eq++)
		{
			m_viscous_flux[idx_eq] = new double* [3];
			m_prim_grad[idx_eq] = new double* [3];
			m_viscous_flux_grad[idx_eq] = new double** [3];
			for (int idx_dim = 0; idx_dim < 3; idx_dim++)
			{
				m_viscous_flux_grad[idx_eq][idx_dim] = new double* [3];
			}
		}
		m_data->GetData("dt", m_dt);
		m_data->GetData("temperture", m_temperture);
		m_data->GetData("non_physical", m_non_physical);
		for (int idx_eq = 0; idx_eq < GetEqNum(); idx_eq++)
		{
			m_data->GetData("primitive_" + std::to_string(idx_eq), m_prim[idx_eq]);
			m_data->GetData("conservative_" + std::to_string(idx_eq), m_cons[idx_eq]);
			m_data->GetData("conservative_old_" + std::to_string(idx_eq), m_cons_old[idx_eq]);
			m_data->GetData("residual_" + std::to_string(idx_eq), m_residual[idx_eq]);
			m_data->GetData("limiter_" + std::to_string(idx_eq), m_limiter[idx_eq]);
			for (int idx_dim = 0; idx_dim < 3; idx_dim++)
			{
				m_data->GetData("viscous_flux_" + std::to_string(idx_eq) + "_" + std::to_string(idx_dim), m_viscous_flux[idx_eq][idx_dim]);
				m_data->GetData("prim_grad_" + std::to_string(idx_eq) + "_" + std::to_string(idx_dim), m_prim_grad[idx_eq][idx_dim]);
				for (int jDim = 0; jDim < 3; jDim++)
				{
					m_data->GetData("viscous_flux_grad_" + std::to_string(idx_eq) + "_" + std::to_string(idx_dim) + "_" + std::to_string(jDim), m_viscous_flux_grad[idx_eq][idx_dim][jDim]);
				}
			}
		}
		m_temperture_grad = new double* [3];
		for (int idx_dim = 0; idx_dim < 3; idx_dim++)
		{
			m_data->GetData("temperture_grad_" + std::to_string(idx_dim), m_temperture_grad[idx_dim]);
		}



	}
}