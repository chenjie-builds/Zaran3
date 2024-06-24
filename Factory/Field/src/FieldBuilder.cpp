#include "FieldBuilder.h"
#include "GridFactory.h"
#include"FNGridFactory.h"
#include"GridStructFactory.h"
#include"NSSolverFNFDM.h"
#include "NSFieldFN.h"
#include "NSFieldStruct.h"
namespace zaran
{
    FieldManager* FieldBuilder::Create()
    {
        FieldManager* global_field = new FieldManager();
        GridFactory* grid_factory;
        if (m_grid_type == GridType::Flexible)
        {
            grid_factory = new FNGridFactorySYSU();
        }
        else if (m_grid_type == GridType::Structured)
        {
            grid_factory = new GridStructFactory();
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
        Field* field;
        if (m_solver_type == FieldSolverType::NS_FNFDM)
            field = new NSFieldFNFDM(grid);
        else if (m_solver_type == FieldSolverType::NS_Struct)
        {
            field = new NSFieldStruct(grid);
        }
        else
        {
            Log::warn("Unsupported Solver Type! Please Check!");
            system("pause");
        }
        global_field->AddField(field);
        delete[] grid_factory;
        return global_field;
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