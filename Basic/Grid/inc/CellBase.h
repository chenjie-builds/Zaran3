#pragma once
namespace zaran
{
class CellBase
{
  public:
    CellBase();
    virtual ~CellBase();
    const int &GetCellNum() const;

  protected:
    void SetCellNum(int cell_num);
  private:
    int m_cell_num;
};

} // namespace zaran