#include "NSSolverStruct.h"
#include "GlobalData.h"
#include "log.h"
#include"MathBasic.h"
namespace zaran
{

	void NSSolverStruct::ComputeCoordTrans()
	{
		auto& grid = GetGrid();
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j, int k) {return grid->GetNodeIndex(i, j, k); };
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeType = nodeTopo->GetType();
		auto& coordTransCoef = m_CoordTrans;
		CoordTrans coordTrans;
		for (int k = ks; k < ke; k++)
		{
			for (int j = js; j < je; j++)
			{
				for (int i = is; i < ie; i++)
				{
					iNode = NodeIndex(i, j, k);
					coordTrans.CalcCoordTrans(3, nodeCoord[NodeIndex(i + 1, j, k)], nodeCoord[NodeIndex(i - 1, j, k)], nodeCoord[NodeIndex(i, j + 1, k)], nodeCoord[NodeIndex(i, j - 1, k)], nodeCoord[NodeIndex(i, j, k + 1)], nodeCoord[NodeIndex(i, j, k - 1)]);
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
	}

	void NSSolverStruct::ComputeGradientWLS()
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
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j, int k) {return grid->GetNodeIndex(i, j, k); };
		Matrix3d A;
		DVector3D b, grad;
		double omega = 0;
		double deltaVal;
		double deltaX, deltaY, deltaZ;
		IArray neighborNodeIndex(6);
		Matrix deltaCoord(3, 6);
		for (int k = ks; k < ke; k++)
		{
			for (int j = js; j < je; j++)
			{
				for (int i = is; i < ie; i++)
				{
					iNode = NodeIndex(i, j, k);
					auto& currentCoord = nodeCoord[iNode];
					neighborNodeIndex[0] = NodeIndex(i + 1, j, k);
					neighborNodeIndex[1] = NodeIndex(i - 1, j, k);
					neighborNodeIndex[2] = NodeIndex(i, j + 1, k);
					neighborNodeIndex[3] = NodeIndex(i, j - 1, k);
					neighborNodeIndex[4] = NodeIndex(i, j, k + 1);
					neighborNodeIndex[5] = NodeIndex(i, j, k - 1);
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
							deltaZ = deltaCoord(2, iNeib);
							A(0, 0) += omega * omega * deltaX * deltaX;
							A(0, 1) += omega * omega * deltaX * deltaY;
							A(0, 2) += omega * deltaX * deltaZ;
							A(1, 0) += omega * omega * deltaY * deltaX;
							A(1, 1) += omega * omega * deltaY * deltaY;
							A(1, 2) += omega * deltaY * deltaZ;
							A(2, 0) += omega * deltaZ * deltaX;
							A(2, 1) += omega * deltaZ * deltaY;
							A(2, 2) += omega * deltaZ * deltaZ;
							b(0) += omega * omega * deltaVal * deltaX;
							b(1) += omega * omega * deltaVal * deltaY;
							b(2) += omega * deltaVal * deltaZ;
						}
						grad = A.inverse() * b;
						(*primGradX[iVal])[iNode] = grad.x();
						(*primGradY[iVal])[iNode] = grad.y();
						(*primGradZ[iVal])[iNode] = grad.z();
					}
				}
			}
		}
	}

	Ptr<StructGrid> NSSolverStruct::NSSolverStruct::GetGrid()
	{
		return std::static_pointer_cast<StructGrid> (Solver::GetGrid());
	}
	void NSSolverStruct::InviscidFlux()
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
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j, int k) {return grid->GetNodeIndex(i, j, k); };
		DVector3D r, grad;
		Ptr<RiemannSolverPara >riemanPara = std::make_shared<RiemannSolverPara>();
		riemanPara->gammaL = riemanPara->gammaR = 1.4;
		for (int k = ks; k < ke; k++)
		{
			for (int j = js; j < je; j++)
			{
				for (int i = is; i < ie; i++)
				{
					iNode = NodeIndex(i, j, k);
					auto& jacobi = (*coordTrans[32])[iNode];
					// i direction
					riemanPara->norm(0) = (*coordTrans[16])[iNode];
					riemanPara->norm(1) = (*coordTrans[17])[iNode];
					riemanPara->norm(2) = (*coordTrans[18])[iNode];
					riemanPara->nt = (*coordTrans[19])[iNode];
					r = nodeCoord[NodeIndex(i + 1, j, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = (*primGradX[iVal])[iNode];
						grad(1) = (*primGradY[iVal])[iNode];
						grad(2) = (*primGradZ[iVal])[iNode];
						riemanPara->primL(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
						grad(0) = (*primGradX[iVal])[NodeIndex(i + 1, j, k)];
						grad(1) = (*primGradY[iVal])[NodeIndex(i + 1, j, k)];
						grad(2) = (*primGradZ[iVal])[NodeIndex(i + 1, j, k)];
						riemanPara->primR(iVal) = (*prim[iVal])[NodeIndex(i + 1, j, k)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemanPara);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
					r = nodeCoord[NodeIndex(i - 1, j, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = (*primGradX[iVal])[iNode];
						grad(1) = (*primGradY[iVal])[iNode];
						grad(2) = (*primGradZ[iVal])[iNode];
						riemanPara->primR(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
						grad(0) = (*primGradX[iVal])[NodeIndex(i - 1, j, k)];
						grad(1) = (*primGradY[iVal])[NodeIndex(i - 1, j, k)];
						grad(2) = (*primGradZ[iVal])[NodeIndex(i - 1, j, k)];
						riemanPara->primL(iVal) = (*prim[iVal])[NodeIndex(i - 1, j, k)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemanPara);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;

					// j direction
					riemanPara->norm(0) = (*coordTrans[20])[iNode];
					riemanPara->norm(1) = (*coordTrans[21])[iNode];
					riemanPara->norm(2) = (*coordTrans[22])[iNode];
					riemanPara->nt = (*coordTrans[23])[iNode];
					r = nodeCoord[NodeIndex(i, j + 1, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = (*primGradX[iVal])[iNode];
						grad(1) = (*primGradY[iVal])[iNode];
						grad(2) = (*primGradZ[iVal])[iNode];
						riemanPara->primL(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
						grad(0) = (*primGradX[iVal])[NodeIndex(i, j + 1, k)];
						grad(1) = (*primGradY[iVal])[NodeIndex(i, j + 1, k)];
						grad(2) = (*primGradZ[iVal])[NodeIndex(i, j + 1, k)];
						riemanPara->primR(iVal) = (*prim[iVal])[NodeIndex(i + 1, j, k)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemanPara);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
					r = nodeCoord[NodeIndex(i, j - 1, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = (*primGradX[iVal])[iNode];
						grad(1) = (*primGradY[iVal])[iNode];
						grad(2) = (*primGradZ[iVal])[iNode];
						riemanPara->primR(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
						grad(0) = (*primGradX[iVal])[NodeIndex(i, j - 1, k)];
						grad(1) = (*primGradY[iVal])[NodeIndex(i, j - 1, k)];
						grad(2) = (*primGradZ[iVal])[NodeIndex(i, j - 1, k)];
						riemanPara->primL(iVal) = (*prim[iVal])[NodeIndex(i - 1, j, k)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemanPara);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;

					// k direction
					riemanPara->norm(0) = (*coordTrans[24])[iNode];
					riemanPara->norm(1) = (*coordTrans[25])[iNode];
					riemanPara->norm(2) = (*coordTrans[26])[iNode];
					riemanPara->nt = (*coordTrans[27])[iNode];
					r = nodeCoord[NodeIndex(i, j, k + 1)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = (*primGradX[iVal])[iNode];
						grad(1) = (*primGradY[iVal])[iNode];
						grad(2) = (*primGradZ[iVal])[iNode];
						riemanPara->primL(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
						grad(0) = (*primGradX[iVal])[NodeIndex(i, j, k + 1)];
						grad(1) = (*primGradY[iVal])[NodeIndex(i, j, k + 1)];
						grad(2) = (*primGradZ[iVal])[NodeIndex(i, j, k + 1)];
						riemanPara->primR(iVal) = (*prim[iVal])[NodeIndex(i + 1, j, k)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemanPara);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
					r = nodeCoord[NodeIndex(i, j, k - 1)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = (*primGradX[iVal])[iNode];
						grad(1) = (*primGradY[iVal])[iNode];
						grad(2) = (*primGradZ[iVal])[iNode];
						riemanPara->primR(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
						grad(0) = (*primGradX[iVal])[NodeIndex(i, j, k - 1)];
						grad(1) = (*primGradY[iVal])[NodeIndex(i, j, k - 1)];
						grad(2) = (*primGradZ[iVal])[NodeIndex(i, j, k - 1)];
						riemanPara->primL(iVal) = (*prim[iVal])[NodeIndex(i - 1, j, k)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemanPara);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;
				}
			}
		}
	}

	void NSSolverStruct::ComputeTimeStepLocal()
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
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j, int k) {return grid->GetNodeIndex(i, j, k); };
		double minDt = LARGE_NUMBER;
		double gamma = 1.4;
		for (int k = ks; k < ke; k++)
		{
			for (int j = js; j < je; j++)
			{
				for (int i = is; i < ie; i++)
				{
					iNode = NodeIndex(i, j, k);
					double c = sqrt(gamma * p[iNode] / rho[iNode]);
					double normXi = sqrt((*coordTrans[16])[iNode] * (*coordTrans[16])[iNode] + (*coordTrans[17])[iNode] * (*coordTrans[17])[iNode] + (*coordTrans[18])[iNode] * (*coordTrans[19])[iNode]);
					double normEta = sqrt((*coordTrans[20])[iNode] * (*coordTrans[20])[iNode] + (*coordTrans[21])[iNode] * (*coordTrans[21])[iNode] + (*coordTrans[22])[iNode] * (*coordTrans[22])[iNode]);
					double normZeta = sqrt((*coordTrans[24])[iNode] * (*coordTrans[24])[iNode] + (*coordTrans[25])[iNode] * (*coordTrans[25])[iNode] + (*coordTrans[26])[iNode] * (*coordTrans[26])[iNode]);
					double uXi = u[iNode] * (*coordTrans[16])[iNode] + v[iNode] * (*coordTrans[17])[iNode] + w[iNode] * (*coordTrans[18])[iNode] + (*coordTrans[19])[iNode];
					double uEta = u[iNode] * (*coordTrans[20])[iNode] + v[iNode] * (*coordTrans[21])[iNode] + w[iNode] * (*coordTrans[22])[iNode] + (*coordTrans[23])[iNode];
					double uZeta = u[iNode] * (*coordTrans[24])[iNode] + v[iNode] * (*coordTrans[25])[iNode] + w[iNode] * (*coordTrans[26])[iNode] + (*coordTrans[27])[iNode];
					double lamda = abs(uXi) + abs(uEta) + abs(uZeta) + c * (normXi + normEta + normZeta);
					dt[iNode] = cfl / lamda;
					minDt = Min(minDt, dt[iNode]);
				}
			}
		}
		GlobalData::Update("dt", minDt);
	}
	void NSSolverStruct::ComputeLimiterCoef()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_Primtive;
		auto& limiterCoef = m_LimiterCoef;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j, int k) {return grid->GetNodeIndex(i, j, k); };
		double maxVal, minVal;
		IArray neighborNodeIndex(6);
		for (int k = ks; k < ke; k++)
		{
			for (int j = js; j < je; j++)
			{
				for (int i = is; i < ie; i++)
				{
					iNode = NodeIndex(i, j, k);
					neighborNodeIndex[0] = NodeIndex(i + 1, j, k);
					neighborNodeIndex[1] = NodeIndex(i - 1, j, k);
					neighborNodeIndex[2] = NodeIndex(i, j + 1, k);
					neighborNodeIndex[3] = NodeIndex(i, j - 1, k);
					neighborNodeIndex[4] = NodeIndex(i, j, k + 1);
					neighborNodeIndex[5] = NodeIndex(i, j, k - 1);
					for (size_t iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						maxVal = minVal = (*prim[iVal])[iNode];
						for (size_t iNeib = 0; iNeib < neighborNodeIndex.size(); ++iNeib)
						{
							maxVal = Max(maxVal, (*prim[iVal])[neighborNodeIndex[iNeib]]);
							minVal = Min(minVal, (*prim[iVal])[neighborNodeIndex[iNeib]]);
						}
						double gradx = (*primGradX[iVal])[iNode];
						double grady = (*primGradY[iVal])[iNode];
						double gradz = (*primGradZ[iVal])[iNode];
						double deltaMax = maxVal - (*prim[iVal])[iNode];
						double deltaMin = minVal - (*prim[iVal])[iNode];
						double tempCoef = LARGE_NUMBER;
						(*limiterCoef[iVal])[iNode] = LARGE_NUMBER;
						for (int iNeighbor = 0; iNeighbor < neighborNodeIndex.size(); ++iNeighbor)
						{
							auto current2Neighbor = nodeCoord[neighborNodeIndex[iNeighbor]] - nodeCoord[iNode];
							double delta2 = 0.5 * (current2Neighbor(0) * gradx + current2Neighbor(1) * grady + current2Neighbor(2) * gradz);
							if (delta2 > 0)
							{
								tempCoef = limiter(maxVal - (*prim[iVal])[iNode], delta2);
							}
							else if (delta2 < 0)
							{
								tempCoef = limiter(minVal - (*prim[iVal])[iNode], delta2);
							}
							else
							{
								tempCoef = 1.0;
							}
							(*limiterCoef[iVal])[iNode] = Min((*limiterCoef[iVal])[iNode], tempCoef);
						}
					}
				}
			}
		}

	}
}