#include "Solver_NS_3D_Struct.h"
#include "GlobalData.h"
#include "log.h"
#include"MathBasic.h"
namespace zaran
{

	void Solver_NS_3D_Struct::CalcMetric()
	{
		Grid_Struct_3D* grid = GetGrid();
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je, ks, ke;
		grid->GetRange(is, ie, js, je, ks, ke);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j, int k) {return grid->GetNodeIndex(i, j, k); };
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeType = nodeTopo->GetType();
		CoordTrans coordTrans;
		for (int k = ks; k < ke; k++)
		{
			for (int j = js; j < je; j++)
			{
				for (int i = is; i < ie; i++)
				{
					iNode = NodeIndex(i, j, k);
					coordTrans.CalcCoordTrans(3, nodeCoord[NodeIndex(i + 1, j, k)], nodeCoord[NodeIndex(i - 1, j, k)], nodeCoord[NodeIndex(i, j + 1, k)], nodeCoord[NodeIndex(i, j - 1, k)], nodeCoord[NodeIndex(i, j, k + 1)], nodeCoord[NodeIndex(i, j, k - 1)]);
					GetMetricXi(iNode)[0] = coordTrans.GetXi()[0];
					GetMetricXi(iNode)[1] = coordTrans.GetXi()[1];
					GetMetricXi(iNode)[2] = coordTrans.GetXi()[2];
					GetMetricXi(iNode)[3] = coordTrans.GetXi()[3];
					GetMetricEta(iNode)[0] = coordTrans.GetEta()[0];
					GetMetricEta(iNode)[1] = coordTrans.GetEta()[1];
					GetMetricEta(iNode)[2] = coordTrans.GetEta()[2];
					GetMetricEta(iNode)[3] = coordTrans.GetEta()[3];
					GetMetricZeta(iNode)[0] = coordTrans.GetZeta()[0];
					GetMetricZeta(iNode)[1] = coordTrans.GetZeta()[1];
					GetMetricZeta(iNode)[2] = coordTrans.GetZeta()[2];
					GetMetricZeta(iNode)[3] = coordTrans.GetZeta()[3];
					GetMetricTau(iNode)[0] = coordTrans.GetTau()[0];
					GetMetricTau(iNode)[1] = coordTrans.GetTau()[1];
					GetMetricTau(iNode)[2] = coordTrans.GetTau()[2];
					GetMetricTau(iNode)[3] = coordTrans.GetTau()[3];
					GetMetricJacob(iNode) = coordTrans.J();
				}
			}
		}
	}

	void Solver_NS_3D_Struct::CalcGradWLS()
	{
		Grid_Struct_3D* grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		auto& prim = m_prim;
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
						GetPrimGrad(iNode, iVal, 0) = grad.x();
						GetPrimGrad(iNode, iVal, 1) = grad.y();
						GetPrimGrad(iNode, iVal, 2) = grad.z();
					}
				}
			}
		}
	}

	Grid_Struct_3D* Solver_NS_3D_Struct::GetGrid()
	{
		return static_cast<Grid_Struct_3D*>(Solver::GetGrid());
	}
	void Solver_NS_3D_Struct::ConvectiveResidual()
	{
		Grid_Struct_3D* grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_prim;
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
					auto& jacobi = GetMetricJacob(iNode);
					// i direction
					riemann_para.norm(0) = GetMetricXi(iNode)[0];
					riemann_para.norm(1) = GetMetricXi(iNode)[1];
					riemann_para.norm(2) = GetMetricXi(iNode)[2];
					riemann_para.nt = GetMetricXi(iNode)[3];
					r = nodeCoord[NodeIndex(i + 1, j, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = GetPrimGrad(iNode, iVal, 0);
						grad(1) = GetPrimGrad(iNode, iVal, 1);
						grad(2) = GetPrimGrad(iNode, iVal, 2);
						riemann_para.prim_left(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
						grad(0) = GetPrimGrad(NodeIndex(i + 1, j, k), iVal, 0);
						grad(1) = GetPrimGrad(NodeIndex(i + 1, j, k), iVal, 1);
						grad(2) = GetPrimGrad(NodeIndex(i + 1, j, k), iVal, 2);
						riemann_para.prim_right(iVal) = m_prim[iVal][NodeIndex(i + 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
					r = nodeCoord[NodeIndex(i - 1, j, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = GetPrimGrad(iNode, iVal, 0);
						grad(1) = GetPrimGrad(iNode, iVal, 1);
						grad(2) = GetPrimGrad(iNode, iVal, 2);
						riemann_para.prim_right(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
						grad(0) = GetPrimGrad(NodeIndex(i - 1, j, k), iVal, 0);
						grad(1) = GetPrimGrad(NodeIndex(i - 1, j, k), iVal, 1);
						grad(2) = GetPrimGrad(NodeIndex(i - 1, j, k), iVal, 2);
						riemann_para.prim_left(iVal) = m_prim[iVal][NodeIndex(i - 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;

					// j direction
					riemann_para.norm(0) = GetMetricEta(iNode)[0];
					riemann_para.norm(1) = GetMetricEta(iNode)[1];
					riemann_para.norm(2) = GetMetricEta(iNode)[2];
					riemann_para.nt = GetMetricEta(iNode)[3];
					r = nodeCoord[NodeIndex(i, j + 1, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = GetPrimGrad(iNode, iVal, 0);
						grad(1) = GetPrimGrad(iNode, iVal, 1);
						grad(2) = GetPrimGrad(iNode, iVal, 2);
						riemann_para.prim_left(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
						grad(0) = GetPrimGrad(NodeIndex(i, j + 1, k), iVal, 0);
						grad(1) = GetPrimGrad(NodeIndex(i, j + 1, k), iVal, 1);
						grad(2) = GetPrimGrad(NodeIndex(i, j + 1, k), iVal, 2);
						riemann_para.prim_right(iVal) = m_prim[iVal][NodeIndex(i + 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
					r = nodeCoord[NodeIndex(i, j - 1, k)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = GetPrimGrad(iNode, iVal, 0);
						grad(1) = GetPrimGrad(iNode, iVal, 1);
						grad(2) = GetPrimGrad(iNode, iVal, 2);
						riemann_para.prim_right(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
						grad(0) = GetPrimGrad(NodeIndex(i, j - 1, k), iVal, 0);
						grad(1) = GetPrimGrad(NodeIndex(i, j - 1, k), iVal, 1);
						grad(2) = GetPrimGrad(NodeIndex(i, j - 1, k), iVal, 2);
						riemann_para.prim_left(iVal) = m_prim[iVal][NodeIndex(i - 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;

					// k direction
					riemann_para.norm(0) = GetMetricZeta(iNode)[0];
					riemann_para.norm(1) = GetMetricZeta(iNode)[1];
					riemann_para.norm(2) = GetMetricZeta(iNode)[2];
					riemann_para.nt = GetMetricZeta(iNode)[3];
					r = nodeCoord[NodeIndex(i, j, k + 1)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = GetPrimGrad(iNode, iVal, 0);
						grad(1) = GetPrimGrad(iNode, iVal, 1);
						grad(2) = GetPrimGrad(iNode, iVal, 2);
						riemann_para.prim_left(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
						grad(0) = GetPrimGrad(NodeIndex(i, j, k + 1), iVal, 0);
						grad(1) = GetPrimGrad(NodeIndex(i, j, k + 1), iVal, 1);
						grad(2) = GetPrimGrad(NodeIndex(i, j, k + 1), iVal, 2);
						riemann_para.prim_right(iVal) = m_prim[iVal][NodeIndex(i + 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
					r = nodeCoord[NodeIndex(i, j, k - 1)] - nodeCoord[NodeIndex(i, j, k)];
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					{
						grad(0) = GetPrimGrad(iNode, iVal, 0);
						grad(1) = GetPrimGrad(iNode, iVal, 1);
						grad(2) = GetPrimGrad(iNode, iVal, 2);
						riemann_para.prim_right(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
						grad(0) = GetPrimGrad(NodeIndex(i, j, k - 1), iVal, 0);
						grad(1) = GetPrimGrad(NodeIndex(i, j, k - 1), iVal, 1);
						grad(2) = GetPrimGrad(NodeIndex(i, j, k - 1), iVal, 2);
						riemann_para.prim_left(iVal) = m_prim[iVal][NodeIndex(i - 1, j, k)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
					}
					riemannSolver_->Solver(riemann_para);
					for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
						GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;
				}
			}
		}
	}

	void Solver_NS_3D_Struct::CalcTimeStepLocal()
	{
		Grid_Struct_3D* grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		FlowSolverPara* para = GetPara();
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
					double normXi = sqrt(GetMetricXi(iNode)[0] * GetMetricXi(iNode)[0] + GetMetricXi(iNode)[1] * GetMetricXi(iNode)[1] + GetMetricXi(iNode)[2] * GetMetricXi(iNode)[2]);
					double normEta = sqrt(GetMetricEta(iNode)[0] * GetMetricEta(iNode)[0] + GetMetricEta(iNode)[1] * GetMetricEta(iNode)[1] + GetMetricEta(iNode)[2] * GetMetricEta(iNode)[2]);
					double normZeta = sqrt(GetMetricZeta(iNode)[0] * GetMetricZeta(iNode)[0] + GetMetricZeta(iNode)[1] * GetMetricZeta(iNode)[1] + GetMetricZeta(iNode)[2] * GetMetricZeta(iNode)[2]);
					double uXi = m_prim[1][iNode] * GetMetricXi(iNode)[0] + m_prim[2][iNode] * GetMetricXi(iNode)[1] + m_prim[3][iNode] * GetMetricXi(iNode)[2];
					double uEta = m_prim[1][iNode] * GetMetricEta(iNode)[0] + m_prim[2][iNode] * GetMetricEta(iNode)[1] + m_prim[3][iNode] * GetMetricEta(iNode)[2];
					double uZeta = m_prim[1][iNode] * GetMetricZeta(iNode)[0] + m_prim[2][iNode] * GetMetricZeta(iNode)[1] + m_prim[3][iNode] * GetMetricZeta(iNode)[2];
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
		Grid_Struct_3D* grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_prim;
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
						double gradx = GetPrimGrad(iNode, iVal, 0);
						double grady = GetPrimGrad(iNode, iVal, 1);
						double gradz = GetPrimGrad(iNode, iVal, 2);
						double deltaMax = maxVal - m_prim[iVal][iNode];
						double deltaMin = minVal - m_prim[iVal][iNode];
						double tempCoef = LARGE_NUMBER;
						GetLimiter(iNode, iVal) = LARGE_NUMBER;
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
							GetLimiter(iNode, iVal) = Min(GetLimiter(iNode, iVal), tempCoef);
						}
					}
				}
			}
		}

	}
}