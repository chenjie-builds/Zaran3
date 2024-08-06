#include "DataManagerNSStruct.h"
namespace zaran
{

    DataManagerNSStruct::DataManagerNSStruct(FieldData* fieldData, int ni, int nj, int nk)
     :DataManagerNS(fieldData, ni* nj* nk)
    {
    }

    DataManagerNSStruct::~DataManagerNSStruct()
    {
        for (int iEqu = 0; iEqu < m_equ_num; iEqu++)
        {
            delete[] m_midnode_prim_left[iEqu];
            delete[] m_midnode_prim_right[iEqu];
        }
        delete[] m_midnode_prim_left;
        delete[] m_midnode_prim_right;
    }

    void DataManagerNSStruct::CreateData()
    {
        DataManagerNS::CreateData();

        FieldDataType type = FieldDataType::real;
        for (int iEqu = 0;iEqu < 5;iEqu++)
        {
            for (int iDim = 0;iDim < 3;iDim++)
            {
                m_data->AddData("midnode_prim_left_" + std::to_string(iEqu) + "_" + std::to_string(iDim), type, m_data_num);
                m_data->AddData("midnode_prim_right_" + std::to_string(iEqu) + "_" + std::to_string(iDim), type, m_data_num);
            }
        }
    }

    void DataManagerNSStruct::RegisterData()
    {
        DataManagerNS::RegisterData();
        m_midnode_prim_left = new double** [m_equ_num];
        m_midnode_prim_right = new double** [m_equ_num];
        for (int iEqu = 0; iEqu < m_equ_num; iEqu++)
        {
            m_midnode_prim_left[iEqu] = new double* [3];
            m_midnode_prim_right[iEqu] = new double* [3];
            for (int iDim = 0; iDim < 3; iDim++)
            {
                m_data->GetData("midnode_prim_left_" + std::to_string(iEqu) + "_" + std::to_string(iDim), m_midnode_prim_left[iEqu][iDim]);
                m_data->GetData("midnode_prim_right_" + std::to_string(iEqu) + "_" + std::to_string(iDim), m_midnode_prim_right[iEqu][iDim]);

            }
        }
    }
    void DataManagerNSStruct::SetMidNodePrimLeft(int iEqu, int iDim, int iNode, double value)
    {
        m_midnode_prim_left[iEqu][iDim][iNode] = value;
    }
    void DataManagerNSStruct::SetMidNodePrimRight(int iEqu, int iDim, int iNode, double value)
    {
        m_midnode_prim_right[iEqu][iDim][iNode] = value;
    }
    void DataManagerNSStruct::SetMidNodePrim(int iEqu, int iDim, int iNode, double value_left, double value_right)
    {
        m_midnode_prim_left[iEqu][iDim][iNode] = value_left;
        m_midnode_prim_right[iEqu][iDim][iNode] = value_right;
    }
    double DataManagerNSStruct::GetMidNodePrimLeft(int iEqu, int iDim, int iNode)
    {
        return m_midnode_prim_left[iEqu][iDim][iNode];
    }
    double DataManagerNSStruct::GetMidNodePrimRight(int iEqu, int iDim, int iNode)
    {
        return m_midnode_prim_right[iEqu][iDim][iNode];
    }
}