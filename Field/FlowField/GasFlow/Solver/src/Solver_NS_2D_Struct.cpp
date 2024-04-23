#include"Solver_NS_2D_Struct.h"
#include"log.h"
#include"MathBasic.h"
namespace zaran
{
	void Solver_NS_2D_Struct::InitData()
	{
		GridPtr grid = GetGrid();
		auto& NodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = NodeTopo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		int initType = para->GetInitFieldType();
		double prim_init[5];
		prim_init[0] = para->GetInflowDensity();
		prim_init[1] = para->GetInflowVelocityX();
		prim_init[2] = para->GetInflowVelocityY();
		prim_init[3] = para->GetInflowVelocityZ();
		prim_init[4] = para->GetInflowPressure();
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double theta = 90 / 180.0 * PI;
		double x, y, z;
		prim_init[0] = 5.4;
		prim_init[1] = 2.2222 * sin(theta);
		prim_init[2] = -2.22222 * cos(theta);
		prim_init[3] = 0;
		prim_init[4] = 10.3333;
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			x = nodeCoord[iNode].x();
			y = nodeCoord[iNode].y();
			if (x >= 0.07 && y <= x * tan(theta) - 0.07 * tan(theta))
			{
				m_prim[0][iNode] = 1.4;
				m_prim[1][iNode] = 0.0;
				m_prim[2][iNode] = 0.0;
				m_prim[3][iNode] = 0.0;
				m_prim[4][iNode] = 1.0;

			}
			else
			{
				m_prim[0][iNode] = 5.4;
				m_prim[1][iNode] = 2.22222 * sin(theta);
				m_prim[2][iNode] = -2.22222 * cos(theta);
				m_prim[3][iNode] = 0.0;
				m_prim[4][iNode] = 10.33333;
			}
		}
		Prim2Cons();
	}
	void Solver_NS_2D_Struct::CalcMetric()
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
		CoordTrans metric;
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				metric.CalcCoordTrans(int(grid->GetDimension()), nodeCoord[NodeIndex(i + 1, j)], nodeCoord[NodeIndex(i - 1, j)], nodeCoord[NodeIndex(i, j + 1)], nodeCoord[NodeIndex(i, j - 1)]);
				for (int iDim = 0; iDim < 4; ++iDim)
				{
					GetMetricXi(iNode)[iDim] = metric.GetXi()[iDim];
					GetMetricEta(iNode)[iDim] = metric.GetEta()[iDim];
					GetMetricZeta(iNode)[iDim] = metric.GetZeta()[iDim];
					GetMetricTau(iNode)[iDim] = metric.GetTau()[iDim];
				}
				GetMetricJacob(iNode) = metric.J();

			}
		}
	}
	void Solver_NS_2D_Struct::CalcGradWLS()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		double cfl = para->GetCflNumber();
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je;
		grid->GetRange(is, ie, js, je);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j) {return grid->GetNodeIndex(i, j); };
		Matrix2d A;
		DVector2D b, grad;
		double omega = 0;
		double deltaVal;
		double deltaX, deltaY, deltaZ;
		IArray neighborNodeIndex(4);
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
				for (size_t iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					A.setZero();
					b.setZero();
					for (size_t iNeib = 0; iNeib < neighborNodeIndex.size(); ++iNeib)
					{
						omega = DistanceOfTwoPoints(nodeCoord[neighborNodeIndex[iNeib]].data(), nodeCoord[iNode].data());
						omega = 1.0 / omega;
						deltaVal = m_prim[iVal][neighborNodeIndex[iNeib]] - m_prim[iVal][iNode];
						deltaX = nodeCoord[neighborNodeIndex[iNeib]].x() - nodeCoord[iNode].x();
						deltaY = nodeCoord[neighborNodeIndex[iNeib]].y() - nodeCoord[iNode].y();
						A(0, 0) += omega * deltaX * deltaX;
						A(0, 1) += omega * deltaX * deltaY;
						A(1, 0) += omega * deltaY * deltaX;
						A(1, 1) += omega * deltaY * deltaY;
						b(0) += omega * deltaVal * deltaX;
						b(1) += omega * deltaVal * deltaY;
					}
					grad = A.inverse() * b;
					GetPrimGrad(iNode, iVal, 0) = grad(0);
					GetPrimGrad(iNode, iVal, 1) = grad(1);
				}
			}
		}
	}
	Ptr<Grid_Struct_2D> Solver_NS_2D_Struct::GetGrid()
	{
		return std::static_pointer_cast<Grid_Struct_2D>(Solver::GetGrid());
	}
	void Solver_NS_2D_Struct::CalcTimeStepLocal()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		FlowSolverPara* para = GetPara();
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

				double c = sqrt(gamma * m_prim[4][iNode] / m_prim[0][iNode]);
				double normXi = sqrt(GetMetricXi(iNode)[0] * GetMetricXi(iNode)[0] + GetMetricXi(iNode)[1] * GetMetricXi(iNode)[1] + GetMetricXi(iNode)[2] * GetMetricXi(iNode)[2]);
				double normEta = sqrt(GetMetricEta(iNode)[0] * GetMetricEta(iNode)[0] + GetMetricEta(iNode)[1] * GetMetricEta(iNode)[1] + GetMetricEta(iNode)[2] * GetMetricEta(iNode)[2]);
				double uXi = m_prim[1][iNode] * GetMetricXi(iNode)[0] + m_prim[2][iNode] * GetMetricXi(iNode)[1] + m_prim[3][iNode] * GetMetricXi(iNode)[2] + GetMetricXi(iNode)[3];
				double uEta = m_prim[1][iNode] * GetMetricEta(iNode)[0] + m_prim[2][iNode] * GetMetricEta(iNode)[1] + m_prim[3][iNode] * GetMetricEta(iNode)[2] + GetMetricEta(iNode)[3];
				double lamda = abs(uXi) + abs(uEta) + c * (normXi + normEta);
				m_dt[iNode] = cfl / lamda;
				minDt = Min(minDt, m_dt[iNode]);
			}
		}
		GlobalData::Update("dt", minDt);
	}
	void Solver_NS_2D_Struct::ConvectiveResidual()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je;
		grid->GetRange(is, ie, js, je);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j) {return grid->GetNodeIndex(i, j); };
		DVector2D r, grad;
		RiemannSolverPara riemann_para;
		riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				auto& jacobi = GetMetricJacob(iNode);
				// i direction
				riemann_para.norm(0) = GetMetricXi(iNode)[0];
				riemann_para.norm(1) = GetMetricXi(iNode)[1];
				riemann_para.norm(2) = GetMetricXi(iNode)[2];
				riemann_para.nt = GetMetricXi(iNode)[3];
				r[0] = nodeCoord[NodeIndex(i + 1, j)][0] - nodeCoord[NodeIndex(i, j)][0];
				r[1] = nodeCoord[NodeIndex(i + 1, j)][1] - nodeCoord[NodeIndex(i, j)][1];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = GetPrimGrad(iNode, iVal, 0);
					grad(1) = GetPrimGrad(iNode, iVal, 1);
					riemann_para.prim_left(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
					grad(0) = GetPrimGrad(NodeIndex(i + 1, j), iVal, 0);
					grad(1) = GetPrimGrad(NodeIndex(i + 1, j), iVal, 1);
					riemann_para.prim_right(iVal) = m_prim[iVal][NodeIndex(i + 1, j)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
				r[0] = nodeCoord[NodeIndex(i - 1, j)][0] - nodeCoord[NodeIndex(i, j)][0];
				r[1] = nodeCoord[NodeIndex(i - 1, j)][1] - nodeCoord[NodeIndex(i, j)][1];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = GetPrimGrad(iNode, iVal, 0);
					grad(1) = GetPrimGrad(iNode, iVal, 1);
					riemann_para.prim_right(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
					grad(0) = GetPrimGrad(NodeIndex(i - 1, j), iVal, 0);
					grad(1) = GetPrimGrad(NodeIndex(i - 1, j), iVal, 1);
					riemann_para.prim_left(iVal) = m_prim[iVal][NodeIndex(i - 1, j)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;

				// j direction
				riemann_para.norm(0) = GetMetricEta(iNode)[0];
				riemann_para.norm(1) = GetMetricEta(iNode)[1];
				riemann_para.norm(2) = GetMetricEta(iNode)[2];
				riemann_para.nt = GetMetricEta(iNode)[3];
				r[0] = nodeCoord[NodeIndex(i, j + 1)][0] - nodeCoord[NodeIndex(i, j)][0];
				r[1] = nodeCoord[NodeIndex(i, j + 1)][1] - nodeCoord[NodeIndex(i, j)][1];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = GetPrimGrad(iNode, iVal, 0);
					grad(1) = GetPrimGrad(iNode, iVal, 1);
					riemann_para.prim_left(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
					grad(0) = GetPrimGrad(NodeIndex(i, j + 1), iVal, 0);
					grad(1) = GetPrimGrad(NodeIndex(i, j + 1), iVal, 1);
					riemann_para.prim_right(iVal) = m_prim[iVal][NodeIndex(i + 1, j)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
				r[0] = nodeCoord[NodeIndex(i, j - 1)][0] - nodeCoord[NodeIndex(i, j)][0];
				r[1] = nodeCoord[NodeIndex(i, j - 1)][1] - nodeCoord[NodeIndex(i, j)][1];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = GetPrimGrad(iNode, iVal, 0);
					grad(1) = GetPrimGrad(iNode, iVal, 1);
					riemann_para.prim_right(iVal) = m_prim[iVal][iNode] /*+ 0.5 * GetLimiter(iNode,iVal) * grad.dot(r)*/;
					grad(0) = GetPrimGrad(NodeIndex(i, j - 1), iVal, 0);
					grad(1) = GetPrimGrad(NodeIndex(i, j - 1), iVal, 1);
					riemann_para.prim_left(iVal) = m_prim[iVal][NodeIndex(i - 1, j)] /*- 0.5 * m_limiter[iVal][NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;
			}
		}
	}
	void Solver_NS_2D_Struct::InviscidFluxMUSCL()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_prim;
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je;
		grid->GetRange(is, ie, js, je);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j) {return grid->GetNodeIndex(i, j); };
		auto Prim = [&](int iVal, int i, int j) {return m_prim[iVal][NodeIndex(i, j)]; };
		RiemannSolverPara riemann_para;
		riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
		double k = GlobalData::GetDouble("muscl_k");
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				auto& jacobi = GetMetricJacob(iNode);
				// i direction
				riemann_para.norm(0) = GetMetricXi(iNode)[0];
				riemann_para.norm(1) = GetMetricXi(iNode)[1];
				riemann_para.norm(2) = GetMetricXi(iNode)[2];
				riemann_para.nt = GetMetricXi(iNode)[3];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.prim_left(iVal) = Prim(iVal, i, j) +
						0.25 * ((1 - k) * (limiter(Prim(iVal, i + 1, j) - Prim(iVal, i, j), Prim(iVal, i, j) - Prim(iVal, i - 1, j)) * (Prim(iVal, i, j) - Prim(iVal, i - 1, j)) +
							(1 + k) * (limiter(Prim(iVal, i, j) - Prim(iVal, i - 1, j), Prim(iVal, i + 1, j) - Prim(iVal, i, j)) * (Prim(iVal, i + 1, j) - Prim(iVal, i, j)))));
					riemann_para.prim_right(iVal) = Prim(iVal, i + 1, j) +
						0.25 * ((1 - k) * (limiter(Prim(iVal, i + 1, j) - Prim(iVal, i, j), Prim(iVal, i + 2, j) - Prim(iVal, i + 1, j)) * (Prim(iVal, i + 2, j) - Prim(iVal, i + 1, j)) +
							(1 + k) * (limiter(Prim(iVal, i + 2, j) - Prim(iVal, i + 1, j), Prim(iVal, i + 1, j) - Prim(iVal, i, j)) * (Prim(iVal, i + 1, j) - Prim(iVal, i, j)))));
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.prim_left(iVal) = Prim(iVal, i - 1, j) +
						0.25 * ((1 - k) * (limiter(Prim(iVal, i, j) - Prim(iVal, i - 1, j), Prim(iVal, i - 1, j) - Prim(iVal, i - 2, j)) * (Prim(iVal, i - 1, j) - Prim(iVal, i - 2, j)) +
							(1 + k) * (limiter(Prim(iVal, i - 1, j) - Prim(iVal, i - 2, j), Prim(iVal, i, j) - Prim(iVal, i - 1, j)) * (Prim(iVal, i, j) - Prim(iVal, i - 1, j)))));
					riemann_para.prim_right(iVal) = Prim(iVal, i, j) +
						0.25 * ((1 - k) * (limiter(Prim(iVal, i, j) - Prim(iVal, i - 1, j), Prim(iVal, i + 1, j) - Prim(iVal, i, j)) * (Prim(iVal, i + 1, j) - Prim(iVal, i, j)) +
							(1 + k) * (limiter(Prim(iVal, i + 1, j) - Prim(iVal, i, j), Prim(iVal, i, j) - Prim(iVal, i - 1, j)) * (Prim(iVal, i, j) - Prim(iVal, i - 1, j)))));
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;

				// j direction
				riemann_para.norm(0) = GetMetricEta(iNode)[0];
				riemann_para.norm(1) = GetMetricEta(iNode)[1];
				riemann_para.norm(2) = GetMetricEta(iNode)[2];
				riemann_para.nt = GetMetricEta(iNode)[3];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.prim_left(iVal) = Prim(iVal, i, j) +
						0.25 * ((1 - k) * (limiter(Prim(iVal, i + 1, j) - Prim(iVal, i, j), Prim(iVal, i, j) - Prim(iVal, i - 1, j)) * (Prim(iVal, i, j) - Prim(iVal, i - 1, j)) +
							(1 + k) * (limiter(Prim(iVal, i, j) - Prim(iVal, i - 1, j), Prim(iVal, i + 1, j) - Prim(iVal, i, j)) * (Prim(iVal, i + 1, j) - Prim(iVal, i, j)))));
					riemann_para.prim_right(iVal) = Prim(iVal, i + 1, j) +
						0.25 * ((1 - k) * (limiter(Prim(iVal, i + 1, j) - Prim(iVal, i, j), Prim(iVal, i + 2, j) - Prim(iVal, i + 1, j)) * (Prim(iVal, i + 2, j) - Prim(iVal, i + 1, j)) +
							(1 + k) * (limiter(Prim(iVal, i + 2, j) - Prim(iVal, i + 1, j), Prim(iVal, i + 1, j) - Prim(iVal, i, j)) * (Prim(iVal, i + 1, j) - Prim(iVal, i, j)))));
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.prim_left(iVal) = Prim(iVal, i - 1, j) +
						0.25 * ((1 - k) * (limiter(Prim(iVal, i, j) - Prim(iVal, i - 1, j), Prim(iVal, i - 1, j) - Prim(iVal, i - 2, j)) * (Prim(iVal, i - 1, j) - Prim(iVal, i - 2, j)) +
							(1 + k) * (limiter(Prim(iVal, i - 1, j) - Prim(iVal, i - 2, j), Prim(iVal, i, j) - Prim(iVal, i - 1, j)) * (Prim(iVal, i, j) - Prim(iVal, i - 1, j)))));
					riemann_para.prim_right(iVal) = Prim(iVal, i, j) +
						0.25 * ((1 - k) * (limiter(Prim(iVal, i, j) - Prim(iVal, i - 1, j), Prim(iVal, i + 1, j) - Prim(iVal, i, j)) * (Prim(iVal, i + 1, j) - Prim(iVal, i, j)) +
							(1 + k) * (limiter(Prim(iVal, i + 1, j) - Prim(iVal, i, j), Prim(iVal, i, j) - Prim(iVal, i - 1, j)) * (Prim(iVal, i, j) - Prim(iVal, i - 1, j)))));
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;
			}
		}
	}
	void Solver_NS_2D_Struct::CalcLimiter()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je;
		grid->GetRange(is, ie, js, je);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j) {return grid->GetNodeIndex(i, j); };
		double maxVal, minVal;
		IArray neighborNodeIndex(4);
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				neighborNodeIndex[0] = NodeIndex(i + 1, j);
				neighborNodeIndex[1] = NodeIndex(i - 1, j);
				neighborNodeIndex[2] = NodeIndex(i, j + 1);
				neighborNodeIndex[3] = NodeIndex(i, j - 1);
				for (size_t iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					maxVal = minVal = m_prim[iVal][iNode];
					for (size_t iNeib = 0; iNeib < neighborNodeIndex.size(); ++iNeib)
					{
						maxVal = Max(maxVal, m_prim[iVal][neighborNodeIndex[iNeib]]);
						minVal = Min(minVal, m_prim[iVal][neighborNodeIndex[iNeib]]);
					}
					double deltaMax = maxVal - m_prim[iVal][iNode];
					double deltaMin = minVal - m_prim[iVal][iNode];
					double tempCoef = LARGE_NUMBER;
					GetLimiter(iNode, iVal) = LARGE_NUMBER;
					for (int iNeighbor = 0; iNeighbor < neighborNodeIndex.size(); ++iNeighbor)
					{
						auto current2Neighbor = nodeCoord[neighborNodeIndex[iNeighbor]] - nodeCoord[iNode];
						double delta2 = 0.5 * (current2Neighbor(0) * GetPrimGrad(iNode, iVal, 0) + current2Neighbor(1) * GetPrimGrad(iNode, iVal)[0]);
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
	void Solver_NS_2D_Struct::InletBC(Boundary& bound)
	{
		FlowSolverPara* para = GetPara();
		int  ghostIndex = bound.GetGhostIndex();
		double prim_init[5];
		prim_init[0] = para->GetInflowDensity();
		prim_init[1] = para->GetInflowVelocityX();
		prim_init[2] = para->GetInflowVelocityY();
		prim_init[3] = para->GetInflowVelocityZ();
		prim_init[4] = para->GetInflowPressure();
		m_prim[0][ghostIndex] = prim_init[0];
		m_prim[1][ghostIndex] = prim_init[1];
		m_prim[2][ghostIndex] = prim_init[2];
		m_prim[3][ghostIndex] = prim_init[3];
		m_prim[4][ghostIndex] = prim_init[4];
		Prim2Cons(m_prim[0][ghostIndex], m_prim[1][ghostIndex], m_prim[2][ghostIndex], m_prim[3][ghostIndex], m_prim[4][ghostIndex], GetCons(0, ghostIndex), GetCons(1, ghostIndex), GetCons(2, ghostIndex), GetCons(3, ghostIndex), GetCons(4, ghostIndex));
	}
	void Solver_NS_2D_Struct::OutletBC(Boundary& bound)
	{
		int boundIndex = bound.GetIndex();
		int ghostIndex = bound.GetGhostIndex();
		m_prim[0][ghostIndex] = m_prim[0][boundIndex];
		m_prim[1][ghostIndex] = m_prim[1][boundIndex];
		m_prim[2][ghostIndex] = m_prim[2][boundIndex];
		m_prim[3][ghostIndex] = m_prim[3][boundIndex];
		m_prim[4][ghostIndex] = m_prim[4][boundIndex];
		Prim2Cons(m_prim[0][ghostIndex], m_prim[1][ghostIndex], m_prim[2][ghostIndex], m_prim[3][ghostIndex], m_prim[4][ghostIndex], GetCons(0, ghostIndex), GetCons(1, ghostIndex), GetCons(2, ghostIndex), GetCons(3, ghostIndex), GetCons(4, ghostIndex));
	}
}