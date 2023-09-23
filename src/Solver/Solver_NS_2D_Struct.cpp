#include"Solver_NS_2D_Struct.h"
#include"log.h"
#include"MathBasic.h"
namespace zaran
{
	void Solver_NS_2D_Struct::ComputeCoordTrans()
	{
		auto& grid = GetGrid();
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je;
		grid->GetRange(is, ie, js, je);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j) {return grid->GetNodeIndex(i, j); };
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeType = nodeTopo->GetType();
		auto& coordTransCoef = m_CoordTrans;
		CoordTrans coordTrans;
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				coordTrans.CalcCoordTrans(int(grid->GetDimension()), nodeCoord[NodeIndex(i + 1, j)], nodeCoord[NodeIndex(i - 1, j)], nodeCoord[NodeIndex(i, j + 1)], nodeCoord[NodeIndex(i, j - 1)]);
				(*coordTransCoef[0])[iNode] = coordTrans.GetX()[0];
				(*coordTransCoef[1])[iNode] = coordTrans.GetX()[1];
				(*coordTransCoef[2])[iNode] = coordTrans.GetX()[2];
				(*coordTransCoef[3])[iNode] = coordTrans.GetX()[3];
				(*coordTransCoef[4])[iNode] = coordTrans.GetY()[0];
				(*coordTransCoef[5])[iNode] = coordTrans.GetY()[1];
				(*coordTransCoef[6])[iNode] = coordTrans.GetY()[2];
				(*coordTransCoef[7])[iNode] = coordTrans.GetY()[3];
				(*coordTransCoef[8])[iNode] = coordTrans.GetZ()[0];
				(*coordTransCoef[9])[iNode] = coordTrans.GetZ()[1];
				(*coordTransCoef[10])[iNode] = coordTrans.GetZ()[2];
				(*coordTransCoef[11])[iNode] = coordTrans.GetZ()[3];
				(*coordTransCoef[12])[iNode] = coordTrans.GetT()[0];
				(*coordTransCoef[13])[iNode] = coordTrans.GetT()[1];
				(*coordTransCoef[14])[iNode] = coordTrans.GetT()[2];
				(*coordTransCoef[15])[iNode] = coordTrans.GetT()[3];
				(*coordTransCoef[16])[iNode] = coordTrans.GetXi()[0];
				(*coordTransCoef[17])[iNode] = coordTrans.GetXi()[1];
				(*coordTransCoef[18])[iNode] = coordTrans.GetXi()[2];
				(*coordTransCoef[19])[iNode] = coordTrans.GetXi()[3];
				(*coordTransCoef[20])[iNode] = coordTrans.GetEta()[0];
				(*coordTransCoef[21])[iNode] = coordTrans.GetEta()[1];
				(*coordTransCoef[22])[iNode] = coordTrans.GetEta()[2];
				(*coordTransCoef[23])[iNode] = coordTrans.GetEta()[3];
				(*coordTransCoef[24])[iNode] = coordTrans.GetZeta()[0];
				(*coordTransCoef[25])[iNode] = coordTrans.GetZeta()[1];
				(*coordTransCoef[26])[iNode] = coordTrans.GetZeta()[2];
				(*coordTransCoef[27])[iNode] = coordTrans.GetZeta()[3];
				(*coordTransCoef[28])[iNode] = coordTrans.GetTau()[0];
				(*coordTransCoef[29])[iNode] = coordTrans.GetTau()[1];
				(*coordTransCoef[30])[iNode] = coordTrans.GetTau()[2];
				(*coordTransCoef[31])[iNode] = coordTrans.GetTau()[3];
				(*coordTransCoef[32])[iNode] = coordTrans.J();
			}
		}
	}
	void Solver_NS_2D_Struct::ComputeGradientWLS()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		FlowSolverParaPtr para = GetPara();
		auto& prim = m_Primtive;
		auto& limiterCoef = m_LimiterCoef;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		double cfl = para->GetCflNumber();
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je;
		grid->GetRange(is, ie, js, je);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j) {return grid->GetNodeIndex(i, j); };
		Matrix3d A;
		DVector3D b, grad;
		double omega = 0;
		double deltaVal;
		double deltaX, deltaY, deltaZ;
		IArray neighborNodeIndex(4);
		Matrix deltaCoord(2, 4);
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				auto& currentCoord = nodeCoord[iNode];
				neighborNodeIndex[0] = NodeIndex(i + 1, j);
				neighborNodeIndex[1] = NodeIndex(i - 1, j);
				neighborNodeIndex[2] = NodeIndex(i, j + 1);
				neighborNodeIndex[3] = NodeIndex(i, j - 1);
				for (int iNeib = 0; iNeib < neighborNodeIndex.size(); iNeib++)
					deltaCoord.col(iNeib) = nodeCoord[neighborNodeIndex[iNeib]] - currentCoord;
				for (size_t iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					A.setZero();
					b.setZero();
					for (size_t iNeib = 0; iNeib < neighborNodeIndex.size(); ++iNeib)
					{
						omega = 1.0 / deltaCoord.col(iNeib).norm();
						deltaVal = (*prim[iVal])[neighborNodeIndex[iNeib]] - (*prim[iVal])[iNode];
						deltaX = deltaCoord(0, iNeib);
						deltaY = deltaCoord(1, iNeib);
						A(0, 0) += omega * deltaX * deltaX;
						A(0, 1) += omega * deltaX * deltaY;
						A(1, 0) += omega * deltaY * deltaX;
						A(1, 1) += omega * deltaY * deltaY;
						b(0) += omega * deltaVal * deltaX;
						b(1) += omega * deltaVal * deltaY;
					}
					grad = A.inverse() * b;
					(*primGradX[iVal])[iNode] = grad.x();
					(*primGradY[iVal])[iNode] = grad.y();
				}
			}
		}
	}
	Ptr<Grid_Struct_2D> Solver_NS_2D_Struct::GetGrid()
	{
		return std::static_pointer_cast<Grid_Struct_2D>(Solver::GetGrid());
	}
	void Solver_NS_2D_Struct::ComputeTimeStepLocal()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		FlowSolverParaPtr para = GetPara();
		auto& rho = *m_Primtive[0];
		auto& u = *m_Primtive[1];
		auto& v = *m_Primtive[2];
		auto& w = *m_Primtive[3];
		auto& p = *m_Primtive[4];
		auto& coordTrans = m_CoordTrans;
		auto& dt = *m_TimeStep;
		double cfl = para->GetCflNumber();
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je;
		grid->GetRange(is, ie, js, je);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j) {return grid->GetNodeIndex(i, j); };
		double minDt = LARGE_NUMBER;
		double gamma = 1.4;
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				double c = sqrt(gamma * p[iNode] / rho[iNode]);
				double normXi = sqrt((*coordTrans[16])[iNode] * (*coordTrans[16])[iNode] + (*coordTrans[17])[iNode] * (*coordTrans[17])[iNode] + (*coordTrans[18])[iNode] * (*coordTrans[19])[iNode]);
				double normEta = sqrt((*coordTrans[20])[iNode] * (*coordTrans[20])[iNode] + (*coordTrans[21])[iNode] * (*coordTrans[21])[iNode] + (*coordTrans[22])[iNode] * (*coordTrans[22])[iNode]);
				double uXi = u[iNode] * (*coordTrans[16])[iNode] + v[iNode] * (*coordTrans[17])[iNode] + w[iNode] * (*coordTrans[18])[iNode] + (*coordTrans[19])[iNode];
				double uEta = u[iNode] * (*coordTrans[20])[iNode] + v[iNode] * (*coordTrans[21])[iNode] + w[iNode] * (*coordTrans[22])[iNode] + (*coordTrans[23])[iNode];
				double lamda = abs(uXi) + abs(uEta) + c * (normXi + normEta);
				dt[iNode] = cfl / lamda;
				minDt = Min(minDt, dt[iNode]);
			}
		}
		GlobalData::Update("dt", minDt);
	}
	void Solver_NS_2D_Struct::InviscidFlux()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_Primtive;
		auto& cons = m_Conservative;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		auto& limiterCoef = m_LimiterCoef;
		auto& res = m_Residual;
		auto& coordTrans = m_CoordTrans;
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je;
		grid->GetRange(is, ie, js, je);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j) {return grid->GetNodeIndex(i, j); };
		DVector3D r, grad;
		Ptr<RiemannSolverPara >riemanPara = std::make_shared<RiemannSolverPara>();
		riemanPara->gammaL = riemanPara->gammaR = 1.4;
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				auto& jacobi = (*coordTrans[32])[iNode];
				// i direction
				riemanPara->norm(0) = (*coordTrans[16])[iNode];
				riemanPara->norm(1) = (*coordTrans[17])[iNode];
				riemanPara->norm(2) = (*coordTrans[18])[iNode];
				riemanPara->nt = (*coordTrans[19])[iNode];
				r = nodeCoord[NodeIndex(i + 1, j)] - nodeCoord[NodeIndex(i, j)];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = (*primGradX[iVal])[iNode];
					grad(1) = (*primGradY[iVal])[iNode];
					grad(2) = (*primGradZ[iVal])[iNode];
					riemanPara->primL(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
					grad(0) = (*primGradX[iVal])[NodeIndex(i + 1, j)];
					grad(1) = (*primGradY[iVal])[NodeIndex(i + 1, j)];
					grad(2) = (*primGradZ[iVal])[NodeIndex(i + 1, j)];
					riemanPara->primR(iVal) = (*prim[iVal])[NodeIndex(i + 1, j)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemanPara);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
				r = nodeCoord[NodeIndex(i - 1, j)] - nodeCoord[NodeIndex(i, j)];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = (*primGradX[iVal])[iNode];
					grad(1) = (*primGradY[iVal])[iNode];
					grad(2) = (*primGradZ[iVal])[iNode];
					riemanPara->primR(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
					grad(0) = (*primGradX[iVal])[NodeIndex(i - 1, j)];
					grad(1) = (*primGradY[iVal])[NodeIndex(i - 1, j)];
					grad(2) = (*primGradZ[iVal])[NodeIndex(i - 1, j)];
					riemanPara->primL(iVal) = (*prim[iVal])[NodeIndex(i - 1, j)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemanPara);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;

				// j direction
				riemanPara->norm(0) = (*coordTrans[20])[iNode];
				riemanPara->norm(1) = (*coordTrans[21])[iNode];
				riemanPara->norm(2) = (*coordTrans[22])[iNode];
				riemanPara->nt = (*coordTrans[23])[iNode];
				r = nodeCoord[NodeIndex(i, j + 1)] - nodeCoord[NodeIndex(i, j)];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = (*primGradX[iVal])[iNode];
					grad(1) = (*primGradY[iVal])[iNode];
					grad(2) = (*primGradZ[iVal])[iNode];
					riemanPara->primL(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
					grad(0) = (*primGradX[iVal])[NodeIndex(i, j + 1)];
					grad(1) = (*primGradY[iVal])[NodeIndex(i, j + 1)];
					grad(2) = (*primGradZ[iVal])[NodeIndex(i, j + 1)];
					riemanPara->primR(iVal) = (*prim[iVal])[NodeIndex(i + 1, j)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemanPara);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
				r = nodeCoord[NodeIndex(i, j - 1)] - nodeCoord[NodeIndex(i, j)];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = (*primGradX[iVal])[iNode];
					grad(1) = (*primGradY[iVal])[iNode];
					grad(2) = (*primGradZ[iVal])[iNode];
					riemanPara->primR(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
					grad(0) = (*primGradX[iVal])[NodeIndex(i, j - 1)];
					grad(1) = (*primGradY[iVal])[NodeIndex(i, j - 1)];
					grad(2) = (*primGradZ[iVal])[NodeIndex(i, j - 1)];
					riemanPara->primL(iVal) = (*prim[iVal])[NodeIndex(i - 1, j)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemanPara);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;
			}
		}
	}
}