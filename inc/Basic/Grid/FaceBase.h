/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FaceBase.h
 * \brief Face Topology Base Class
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
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