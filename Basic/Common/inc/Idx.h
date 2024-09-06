#pragma once
class Idx
{
public:
    Idx(size_t idx) {};
    virtual ~Idx() = default;
    virtual void SetIdx(size_t idx) = 0;
    virtual size_t GetIdx() = 0;
};
class Idx1D : public Idx
{
public:
    Idx1D(size_t idx) : Idx(idx) {}
    virtual ~Idx1D() = default;
    void SetIdx(size_t idx) override
    {
        m_idx = idx;
    }
    size_t GetIdx() override
    {
        return m_idx;
    }
private:
    size_t m_idx;
};
class Idx2D : public Idx
{
public:
    Idx2D(size_t idx1, size_t idx2) : Idx(idx1), m_idx2(idx2) {}
    virtual ~Idx2D() = default;
    void SetIdx(size_t idx1, size_t idx2)
    {
        m_idx = idx1;
        m_idx2 = idx2;
    }
    size_t GetIdx1()
    {
        return m_idx;
    }
    size_t GetIdx2()
    {
        return m_idx2;
    }
private:
    size_t m_idx;
    size_t m_idx2;
};