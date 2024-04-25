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
        GridListFactory* grid_factory;
        if (m_dim == Dimension::two)
        {
            if (m_grid_type == GridType::Structured)
            {
                grid_factory = new GridFactoryStruct2D();
            }
            else if (m_grid_type == GridType::Flexible)
            {
                grid_factory = new GridFactoryFNFDM2D();
            }

        }
        else if (m_dim == Dimension::three)
        {
            if (m_grid_type == GridType::Flexible)
            {
                grid_factory = new GridFactoryFNFDM3D();
            }
            else if (m_grid_type == GridType::Zaran)
            {
                grid_factory = new GridFactoryZaran3D();
            }
        }
        else
        {
            Log::warn("Unsupported Dimension! Please Check!");
            system("pause");
        }
        Grid* grid = nullptr;
        grid_factory->Create(grid);
        Solver* solver=nullptr;
        SolverFactory* solver_factory = new SolverFactory();
        solver_factory->Create(grid, solver, m_solver_type);
        m_field = new Field * [1];
        for (int i = 0; i < 1; ++i)
        {
            m_field[i] = new Field();
            m_field[i]->SetGrid(grid);
            m_field[i]->SetSolver(static_cast<FieldSolver*> (solver));
        }
        delete[] grid_factory;
        delete[] solver_factory;
    }

    void FieldFactory::CreateGrid()
    {

    }

    void FieldFactory::CreateField()
    {

    }

    void FieldFactory::CreateSolver()
    {

    }



}