#include"Simulation.h"
#include "log.h"
#include"BasicType.h"

int main()
{
	using namespace zaran;
	Log::Satrt();
	Ptr< Simulation > simulation = std::make_shared<Simulation>();
	simulation->Start();
	ZaranLog::info("计算结束！按任意键退出……");
	system("pause");
}