#include"Test.h"
#include"RiemannSolverPara.h"
#include"Log.h"
namespace zaran
{
    void TestVanLeer()
    {
        RiemannSolverPara para1;
        para1.prim_left << 1.2, 3.0, 0.0, 0.0, 1.0;
        para1.prim_right << 1.2, 3.0, 0.0, 0.0, 1.0;
        para1.gamma_left = 1.4;
        para1.gamma_right = 1.4;
        para1.norm << 0.0, 1.0, 0.0;
        para1.nt = 1.0;
        Vanleer vanleer;
        vanleer.Solver(para1);
        RiemannSolverPara para2;
        para2.prim_left << 1.0, 3.0, 0.0, 0.0, 1.0;
        para2.prim_right << 1.0, 3.0, 0.0, 0.0, 1.0;
        para2.gamma_left = 1.4;
        para2.gamma_right = 1.4;
        para2.norm << 0.0, 1.0, 0.0;
        para2.nt = 1.0;
        vanleer.Solver(para2);
        Log::info("res={},{},{},{},{}", para1.flux(0)-para2.flux(0), para1.flux(1)-para2.flux(1), para1.flux(2)-para2.flux(2), para1.flux(3)-para2.flux(3), para1.flux(4)-para2.flux(4));
    }
}
