#include "NSFieldBlockExplosion.h"
namespace zaran
{


	NSFieldBlockExplosion::NSFieldBlockExplosion(std::shared_ptr<GridBlock>grid) :NSFieldStruct(grid)
	{
	}
	void NSFieldBlockExplosion::AllocateSolver()
	{
		m_solver = std::make_shared<NSSolverBlock>(GetIdx(), "NS_BLOCK_EXPLOSION", GetSolverPara(), GetGrid(), GetDataManager());
	}

	NSFieldBlockExplosion::~NSFieldBlockExplosion()
	{

	}

	void NSFieldBlockExplosion::Allocate()
	{
		NSFieldStruct::Allocate();
	}

	std::shared_ptr<zaran::GridBlock> NSFieldBlockExplosion::GetGrid()
	{
		return std::static_pointer_cast<GridBlock>(NSFieldStruct::GetGrid());
	}

	std::shared_ptr<zaran::NSSolverBlock> NSFieldBlockExplosion::GetSolver()
	{
		return std::static_pointer_cast<NSSolverBlock>(NSFieldStruct::GetSolver());
	}

} // namespace zaran