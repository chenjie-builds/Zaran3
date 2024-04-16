#include "FieldFactory.h"
#include "GridListFactory.h"
#include"Read_Grid_Struct_2D.h"
// #include"Read_Grid_Struct_3D.h"
#include"FNFDM2D.h"
#include"FNFDM3D.h"
#include "Create_Zaran_3D.h"
#include"SolverFactory.h"
namespace zaran
{
    void FieldFactory::Create()
    {
        Ptr<GridListFactory> grid_list_factory;
        Ptr<GridList> grid_list;
        if (m_grid_type == GridType::Structured_2D)
        {
            grid_list_factory = std::make_shared<GridListFactoryStruct2D>();
        }
        // else if(m_grid_type == GridType::Structured_3D)
        // {
        //     grid_list_factory=std::make_shared<GridListFactoryStruct3D>();
        // }
        else if (m_grid_type == GridType::Flexible_2D)
        {
            grid_list_factory = std::make_shared<GridListFactoryFNFDM2D>();
        }
        else if (m_grid_type == GridType::Flexible_3D)
        {
            grid_list_factory = std::make_shared<GridListFactoryFNFDM3D>();
        }
        else if (m_grid_type == GridType::Zaran_3D)
        {
            grid_list_factory = std::make_shared<GridListFactoryZaran3D>();
        }
        else
        {
            Log::warn("Unsupported Grid Type! Please Check!");
            system("pause");
        }
        grid_list_factory->Create(grid_list);
        Ptr<SolverVec> solver_vec;
        Ptr<SolverFactory> solver_factory = std::make_shared<SolverFactory>();
        solver_factory->Create(grid_list, solver_vec, m_solver_type);
        m_field_array.resize(solver_vec->GetSolverNumber());
        for (int i = 0; i < solver_vec->GetSolverNumber(); ++i)
        {
            m_field_array[i] = std::make_shared<Field>();
            m_field_array[i]->SetGrid(grid_list->GetGrid(i));
            m_field_array[i]->SetSolver(std::dynamic_pointer_cast<FieldSolver> (solver_vec->GetSolver(i)));
        }
    }



}