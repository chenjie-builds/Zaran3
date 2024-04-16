#include"Simulation.h"
#include "log.h"
#include"BasicType.h"
#include"Test.h"
int main()
{
	using namespace zaran;
	Logger::Start();
	Ptr< Simulation > simulation = std::make_shared<Simulation>();
	simulation->Run();
	Log::info("计算结束！按任意键退出……");
	system("pause");
}