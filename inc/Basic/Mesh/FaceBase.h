#pragma once
namespace zaran
{
class FaceBase
{
  public:
    FaceBase();
    virtual ~FaceBase();
    const int &GetFaceNum() const;

  protected:
    void SetFaceNum(int face_num);
  private:
    int m_face_num;
};

}