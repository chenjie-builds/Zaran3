#include"Simulation.h"
#include "log.h"
#include"BasicType.h"
#include"Test.h"
int main()
{
	using namespace zaran;
	Logger::Start();
	Simulation* simulation = new Simulation();
	simulation->Run();
	delete[] simulation;
	Log::info("计算结束！按任意键退出……");
	system("pause");
}