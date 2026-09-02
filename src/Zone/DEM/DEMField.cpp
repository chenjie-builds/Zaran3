#include "DEMField.h"

namespace zaran
{

DEMField::DEMField()
    : Field(nullptr, FieldType::DEM)
{
}

void DEMField::Allocate()
{
    AllocateFieldData();
    AllocateSolverPara();
    AllocateSolver();
    AllocateDataManager();
}

void DEMField::AllocateFieldData()
{
    m_dem_data = make_shared<DEMFieldData>();
}

void DEMField::AllocateSolverPara()
{
    m_dem_param = make_shared<DEMSolverParam>();
    m_solver_para = m_dem_param;
}

void DEMField::AllocateSolver()
{
    m_dem_solver = make_shared<DEMSolver>(m_idx, "DEMSolver", m_dem_param, m_dem_data);
    m_solver = m_dem_solver;
}

void DEMField::AllocateDataManager()
{
    m_dem_data_manager = make_shared<DEMDataManager>(m_dem_data);
}

} // namespace zaran
