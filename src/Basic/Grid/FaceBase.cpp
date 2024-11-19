#include "FaceBase.h"
using namespace zaran;
FaceBase::FaceBase()
{
}
FaceBase::~FaceBase()
{
}
const int &FaceBase::GetFaceNum() const
{
    return m_face_num;
}
void FaceBase::SetFaceNum(int face_num)
{
    m_face_num = face_num;
}
