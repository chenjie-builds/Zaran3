#include "FieldBuilder.h"
#include "GridBuilder.h"
#include"FNGridBuilder.h"
#include"NSSolverFNFDM.h"
#include "FieldNS.h"
namespace zaran
{
    void FieldBuilder::Create()
    {
        GridCreater* grid_factory;

        if (m_grid_type == GridType::Flexible)
        {
            grid_factory = new FNGridBuilder();
        }
        // else if (m_grid_type == GridType::Zaran)
        // {
        //     grid_factory = new GridFactoryZaran3D();
        // }

        else
        {
            Log::warn("Unsupported Dimension! Please Check!");
            system("pause");
        }
        GridBase* grid = grid_factory->CreateGrid();
        m_field = new Field * [1];
        for (int i = 0; i < 1; ++i)
        {
            if (m_solver_type == FieldSolverType::NS_FNFDM)
                m_field[i] = new FieldNS_FNFDM(grid);
            else
            {
                Log::warn("Unsupported Solver Type! Please Check!");
                system("pause");
            }
        }
        delete[] grid_factory;
    }

    void FieldBuilder::CreateGrid()
    {

    }

    void FieldBuilder::CreateField()
    {

    }

    void FieldBuilder::CreateSolver()
    {

    }



}