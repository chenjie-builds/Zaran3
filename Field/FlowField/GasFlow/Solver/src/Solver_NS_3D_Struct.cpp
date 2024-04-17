#include "Solver_NS_3D_Struct.h"
#include "GlobalData.h"
#include "log.h"
#include"MathBasic.h"
namespace zaran
{

	void Solver_NS_3D_Struct::CalcMetric()
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
		auto& coordTransCoef = m_metric;
		CoordTrans coordTrans;
		for (int k = ks; k < ke; k++)
		{
			for (int j = js; j < je; j++)
			{
				for (int i = is; i < ie; i++)
				{
					iNode = NodeIndex(i, j, k);
					coordTrans.CalcCoordTrans(3, nodeCoord[NodeIndex(i + 1, j, k)], nodeCoord[NodeIndex(i - 1, j, k)], nodeCoord[NodeIndex(i, j + 1, k)], nodeCoord[NodeIndex(i, j - 1, k)], nodeCoord[NodeIndex(i, j, k + 1)], nodeCoord[NodeIndex(i, j, k - 1)]);
					m_metric[0][iNode] = coordTrans.GetX()[0];
					m_metric[1][iNode] = coordTrans.GetX()[1];
					m_metric[2][iNode] = coordTrans.GetX()[2];
					m_metric[3][iNode] = coordTrans.GetX()[3];
					m_metric[4][iNode] = coordTrans.GetY()[0];
					m_metric[5][iNode] = coordTrans.GetY()[1];
					m_metric[6][iNode] = coordTrans.GetY()[2];
					m_metric[7][iNode] = coordTrans.GetY()[3];
					m_metric[8][iNode] = coordTrans.GetZ()[0];
					m_metric[9][iNode] = coordTrans.GetZ()[1];
					m_metric[10][iNode] = coordTrans.GetZ()[2];
					m_metric[11][iNode] = coordTrans.GetZ()[3];
					m_metric[12][iNode] = coordTrans.GetT()[0];
					m_metric[13][iNode] = coordTrans.GetT()[1];
					m_metric[14][iNode] = coordTrans.GetT()[2];
					m_metric[15][iNode] = coordTrans.GetT()[3];
					m_metric[16][iNode] = coordTrans.GetXi()[0];
					m_metric[17][iNode] = coordTrans.GetXi()[1];
					m_metric[18][iNode] = coordTrans.GetXi()[2];
					m_metric[19][iNode] = coordTrans.GetXi()[3];
					m_metric[20][iNode] = coordTrans.GetEta()[0];
					m_metric[21][iNode] = coordTrans.GetEta()[1];
					m_metric[22][iNode] = coordTrans.GetEta()[2];
					m_metric[23][iNode] = coordTrans.GetEta()[3];
					m_metric[24][iNode] = coordTrans.GetZeta()[0];
					m_metric[25][iNode] = coordTrans.GetZeta()[1];
					m_metric[26][iNode] = coordTrans.GetZeta()[2];
					m_metric[27][iNode] = coordTrans.GetZeta()[3];
					m_metric[28][iNode] = coordTrans.GetTau()[0];
					m_metric[29][iNode] = coordTrans.GetTau()[1];
					m_metric[30][iNode] = coordTrans.GetTau()[2];
					m_metric[31][iNode] = coordTrans.GetTau()[3];
					m_metric[32][iNode] = coordTrans.J();
				}
			}
		}
	}

	void Solver_NS_3D_Struct::CalcGradWLS()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		auto& prim = m_prim;
		auto& limiterCoef = m_limiter;
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
							deltaVal = m_prim[iVal][neighborNodeIndex[iNeib]] - m_prim[iVal][iNode];
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
						m_PrimGradX[iVal][iNode] = grad.x();
						m_PrimGradY[iVal][iNode] = grad.y();
						m_PrimGradZ[iVal][iNode] = grad.z();
					}
				}
			}
		}
	}

	Ptr<Grid_Struct_3D> Solver_NS_3D_Struct::Solver_NS_3D_Struct::GetGrid()
	{
		return std::static_pointer_cast<Grid_Struct_3D> (Solver::GetGrid());
	}
	void Solver_NS_3D_Struct::InviscidFlux()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_prim;
		auto& cons = m_cons;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		auto& limiterCoef = m_limiter;
		auto& res = m_residual;
		auto& coordTrans = m_metric;
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j, int k) {return grid->GetNodeIndex(i, j, k); };
		DVector3D r, grad;
		RiemannSolverPara riemann_para;
		riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
		for (int k = ks; k < ke; k++)
		{
			for (int j = js; j < je; j++)
			{
				for (int i = is; i < ie; i++)
				{
					iNode = NodeIndex(i, j, k);
					auto& jacobi = m_metric[32][iNode];
					// i direction
					riemann_para.norm(0) = m_metric[16][iNode];
					riemann_para.norm(1) = m_metric[17][iNode];
					riemann_para.norm(2) = m_metric[18][iNode];
					riemann_para.nt = m_metric[19][iNode];
					r = nodeCoord[NodeIndex(i + 1, j, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = m_PrimGradX[iVal][iNode];
						grad(1) = m_PrimGradY[iVal][iNode];
						grad(2) = m_PrimGradZ[iVal][iNode];
						riemann_para.prim_left(iVal) = m_prim[iVal][iNode] /*+ 0.5 * m_limiter[iVal][iNode] * grad.dot(r)*/;
						grad(0) = m_PrimGradX[iVal][NodeIndex(i + 1, j, k)];
						grad(1) = m_PrimGradY[iVal][NodeIndex(i + 1, j, k)];
						grad(2) = m_PrimGradZ[iVal][NodeIndex(i + 1, j, k)];
						riemann_para.prim_right(iVal) = m_prim[iVal][NodeIndex(i + 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						m_residual[iVal][iNode] += riemann_para.flux[iVal] / jacobi;
					r = nodeCoord[NodeIndex(i - 1, j, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = m_PrimGradX[iVal][iNode];
						grad(1) = m_PrimGradY[iVal][iNode];
						grad(2) = m_PrimGradZ[iVal][iNode];
						riemann_para.prim_right(iVal) = m_prim[iVal][iNode] /*+ 0.5 * m_limiter[iVal][iNode] * grad.dot(r)*/;
						grad(0) = m_PrimGradX[iVal][NodeIndex(i - 1, j, k)];
						grad(1) = m_PrimGradY[iVal][NodeIndex(i - 1, j, k)];
						grad(2) = m_PrimGradZ[iVal][NodeIndex(i - 1, j, k)];
						riemann_para.prim_left(iVal) = m_prim[iVal][NodeIndex(i - 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						m_residual[iVal][iNode] -= riemann_para.flux[iVal] / jacobi;

					// j direction
					riemann_para.norm(0) = m_metric[20][iNode];
					riemann_para.norm(1) = m_metric[21][iNode];
					riemann_para.norm(2) = m_metric[22][iNode];
					riemann_para.nt = m_metric[23][iNode];
					r = nodeCoord[NodeIndex(i, j + 1, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = m_PrimGradX[iVal][iNode];
						grad(1) = m_PrimGradY[iVal][iNode];
						grad(2) = m_PrimGradZ[iVal][iNode];
						riemann_para.prim_left(iVal) = m_prim[iVal][iNode] /*+ 0.5 * m_limiter[iVal][iNode] * grad.dot(r)*/;
						grad(0) = m_PrimGradX[iVal][NodeIndex(i, j + 1, k)];
						grad(1) = m_PrimGradY[iVal][NodeIndex(i, j + 1, k)];
						grad(2) = m_PrimGradZ[iVal][NodeIndex(i, j + 1, k)];
						riemann_para.prim_right(iVal) = m_prim[iVal][NodeIndex(i + 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						m_residual[iVal][iNode] += riemann_para.flux[iVal] / jacobi;
					r = nodeCoord[NodeIndex(i, j - 1, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = m_PrimGradX[iVal][iNode];
						grad(1) = m_PrimGradY[iVal][iNode];
						grad(2) = m_PrimGradZ[iVal][iNode];
						riemann_para.prim_right(iVal) = m_prim[iVal][iNode] /*+ 0.5 * m_limiter[iVal][iNode] * grad.dot(r)*/;
						grad(0) = m_PrimGradX[iVal][NodeIndex(i, j - 1, k)];
						grad(1) = m_PrimGradY[iVal][NodeIndex(i, j - 1, k)];
						grad(2) = m_PrimGradZ[iVal][NodeIndex(i, j - 1, k)];
						riemann_para.prim_left(iVal) = m_prim[iVal][NodeIndex(i - 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						m_residual[iVal][iNode] -= riemann_para.flux[iVal] / jacobi;

					// k direction
					riemann_para.norm(0) = m_metric[24][iNode];
					riemann_para.norm(1) = m_metric[25][iNode];
					riemann_para.norm(2) = m_metric[26][iNode];
					riemann_para.nt = m_metric[27][iNode];
					r = nodeCoord[NodeIndex(i, j, k + 1)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = m_PrimGradX[iVal][iNode];
						grad(1) = m_PrimGradY[iVal][iNode];
						grad(2) = m_PrimGradZ[iVal][iNode];
						riemann_para.prim_left(iVal) = m_prim[iVal][iNode] /*+ 0.5 * m_limiter[iVal][iNode] * grad.dot(r)*/;
						grad(0) = m_PrimGradX[iVal][NodeIndex(i, j, k + 1)];
						grad(1) = m_PrimGradY[iVal][NodeIndex(i, j, k + 1)];
						grad(2) = m_PrimGradZ[iVal][NodeIndex(i, j, k + 1)];
						riemann_para.prim_right(iVal) = m_prim[iVal][NodeIndex(i + 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						m_residual[iVal][iNode] += riemann_para.flux[iVal] / jacobi;
					r = nodeCoord[NodeIndex(i, j, k - 1)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = m_PrimGradX[iVal][iNode];
						grad(1) = m_PrimGradY[iVal][iNode];
						grad(2) = m_PrimGradZ[iVal][iNode];
						riemann_para.prim_right(iVal) = m_prim[iVal][iNode] /*+ 0.5 * m_limiter[iVal][iNode] * grad.dot(r)*/;
						grad(0) = m_PrimGradX[iVal][NodeIndex(i, j, k - 1)];
						grad(1) = m_PrimGradY[iVal][NodeIndex(i, j, k - 1)];
						grad(2) = m_PrimGradZ[iVal][NodeIndex(i, j, k - 1)];
						riemann_para.prim_left(iVal) = m_prim[iVal][NodeIndex(i - 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						m_residual[iVal][iNode] -= riemann_para.flux[iVal] / jacobi;
				}
			}
		}
	}

	void Solver_NS_3D_Struct::CalcTimeStepLocal()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		FlowSolverPara* para = GetPara();
		auto& coordTrans = m_metric;
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
					double c = sqrt(gamma * m_prim[3][iNode] / m_prim[0][iNode]);
					double normXi = sqrt(m_metric[16][iNode] * m_metric[16][iNode] + m_metric[17][iNode] * m_metric[17][iNode] + m_metric[18][iNode] * m_metric[18][iNode]);
					double normEta = sqrt(m_metric[20][iNode] * m_metric[20][iNode] + m_metric[21][iNode] * m_metric[21][iNode] + m_metric[22][iNode] * m_metric[22][iNode]);
					double normZeta = sqrt(m_metric[24][iNode] * m_metric[24][iNode] + m_metric[25][iNode] * m_metric[25][iNode] + m_metric[26][iNode] * m_metric[26][iNode]);
					double uXi = m_prim[1][iNode] * m_metric[16][iNode] + m_prim[2][iNode] * m_metric[17][iNode] + m_prim[3][iNode] * m_metric[18][iNode];
					double uEta = m_prim[1][iNode] * m_metric[20][iNode] + m_prim[2][iNode] * m_metric[21][iNode] + m_prim[3][iNode] * m_metric[22][iNode];
					double uZeta = m_prim[1][iNode] * m_metric[24][iNode] + m_prim[2][iNode] * m_metric[25][iNode] + m_prim[3][iNode] * m_metric[26][iNode];
					double lambda = abs(uXi) + abs(uEta) + abs(uZeta) + c * (normXi + normEta + normZeta);
					m_dt[iNode] = cfl / lambda;
					minDt = Min(minDt, m_dt[iNode]);
				}
			}
		}
		GlobalData::Update("dt", minDt);
	}
	void Solver_NS_3D_Struct::CalcLimiter()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_prim;
		auto& limiterCoef = m_limiter;
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
						maxVal = minVal = m_prim[iVal][iNode];
						for (size_t iNeib = 0; iNeib < neighborNodeIndex.size(); ++iNeib)
						{
							maxVal = Max(maxVal, m_prim[iVal][neighborNodeIndex[iNeib]]);
							minVal = Min(minVal, m_prim[iVal][neighborNodeIndex[iNeib]]);
						}
						double gradx = m_PrimGradX[iVal][iNode];
						double grady = m_PrimGradY[iVal][iNode];
						double gradz = m_PrimGradZ[iVal][iNode];
						double deltaMax = maxVal - m_prim[iVal][iNode];
						double deltaMin = minVal - m_prim[iVal][iNode];
						double tempCoef = LARGE_NUMBER;
						m_limiter[iVal][iNode] = LARGE_NUMBER;
						for (int iNeighbor = 0; iNeighbor < neighborNodeIndex.size(); ++iNeighbor)
						{
							auto current2Neighbor = nodeCoord[neighborNodeIndex[iNeighbor]] - nodeCoord[iNode];
							double delta2 = 0.5 * (current2Neighbor(0) * gradx + current2Neighbor(1) * grady + current2Neighbor(2) * gradz);
							if (delta2 > 0)
							{
								tempCoef = limiter(maxVal - m_prim[iVal][iNode], delta2);
							}
							else if (delta2 < 0)
							{
								tempCoef = limiter(minVal - m_prim[iVal][iNode], delta2);
							}
							else
							{
								tempCoef = 1.0;
							}
							m_limiter[iVal][iNode] = Min(m_limiter[iVal][iNode], tempCoef);
						}
					}
				}
			}
		}

	}
}