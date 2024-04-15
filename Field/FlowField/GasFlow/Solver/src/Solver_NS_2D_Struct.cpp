#include"Solver_NS_2D_Struct.h"
#include"log.h"
#include"MathBasic.h"
namespace zaran
{
	void Solver_NS_2D_Struct::InitField()
	{
		GridPtr grid = GetGrid();
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
		auto& NodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = NodeTopo->GetCoordinate();
		FlowSolverParaPtr para = GetPara();
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
				rho[iNode] = 1.4;
				u[iNode] = 0.0;
				v[iNode] = 0.0;
				w[iNode] = 0.0;
				p[iNode] = 1.0;
			}
			else
			{
				rho[iNode] = 5.4;
				u[iNode] = 2.22222 * sin(theta);
				v[iNode] = -2.22222 * cos(theta);
				w[iNode] = 0.0;
				p[iNode] = 10.33333;
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
	void Solver_NS_2D_Struct::CalcGradWLS()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		FlowSolverParaPtr para = GetPara();
		auto& prim = m_Primitive;
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
						deltaVal = (*prim[iVal])[neighborNodeIndex[iNeib]] - (*prim[iVal])[iNode];
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
	void Solver_NS_2D_Struct::CalcTimeStepLocal()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		FlowSolverParaPtr para = GetPara();
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
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
		auto& prim = m_Primitive;
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
		DVector2D r, grad;
		RiemannSolverPara riemann_para;
		riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				auto& jacobi = (*coordTrans[32])[iNode];
				// i direction
				riemann_para.norm(0) = (*coordTrans[16])[iNode];
				riemann_para.norm(1) = (*coordTrans[17])[iNode];
				riemann_para.norm(2) = (*coordTrans[18])[iNode];
				riemann_para.nt = (*coordTrans[19])[iNode];
				r[0] = nodeCoord[NodeIndex(i + 1, j)][0] - nodeCoord[NodeIndex(i, j)][0];
				r[1] = nodeCoord[NodeIndex(i + 1, j)][1] - nodeCoord[NodeIndex(i, j)][1];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = (*primGradX[iVal])[iNode];
					grad(1) = (*primGradY[iVal])[iNode];
					riemann_para.prim_left(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
					grad(0) = (*primGradX[iVal])[NodeIndex(i + 1, j)];
					grad(1) = (*primGradY[iVal])[NodeIndex(i + 1, j)];
					riemann_para.prim_right(iVal) = (*prim[iVal])[NodeIndex(i + 1, j)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
				r[0] = nodeCoord[NodeIndex(i - 1, j)][0] - nodeCoord[NodeIndex(i, j)][0];
				r[1] = nodeCoord[NodeIndex(i - 1, j)][1] - nodeCoord[NodeIndex(i, j)][1];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = (*primGradX[iVal])[iNode];
					grad(1) = (*primGradY[iVal])[iNode];
					riemann_para.prim_right(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
					grad(0) = (*primGradX[iVal])[NodeIndex(i - 1, j)];
					grad(1) = (*primGradY[iVal])[NodeIndex(i - 1, j)];
					riemann_para.prim_left(iVal) = (*prim[iVal])[NodeIndex(i - 1, j)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;

				// j direction
				riemann_para.norm(0) = (*coordTrans[20])[iNode];
				riemann_para.norm(1) = (*coordTrans[21])[iNode];
				riemann_para.norm(2) = (*coordTrans[22])[iNode];
				riemann_para.nt = (*coordTrans[23])[iNode];
				r[0] = nodeCoord[NodeIndex(i, j + 1)][0] - nodeCoord[NodeIndex(i, j)][0];
				r[1] = nodeCoord[NodeIndex(i, j + 1)][1] - nodeCoord[NodeIndex(i, j)][1];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = (*primGradX[iVal])[iNode];
					grad(1) = (*primGradY[iVal])[iNode];
					riemann_para.prim_left(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
					grad(0) = (*primGradX[iVal])[NodeIndex(i, j + 1)];
					grad(1) = (*primGradY[iVal])[NodeIndex(i, j + 1)];
					riemann_para.prim_right(iVal) = (*prim[iVal])[NodeIndex(i + 1, j)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
				r[0] = nodeCoord[NodeIndex(i, j - 1)][0] - nodeCoord[NodeIndex(i, j)][0];
				r[1] = nodeCoord[NodeIndex(i, j - 1)][1] - nodeCoord[NodeIndex(i, j)][1];
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					grad(0) = (*primGradX[iVal])[iNode];
					grad(1) = (*primGradY[iVal])[iNode];
					riemann_para.prim_right(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
					grad(0) = (*primGradX[iVal])[NodeIndex(i, j - 1)];
					grad(1) = (*primGradY[iVal])[NodeIndex(i, j - 1)];
					riemann_para.prim_left(iVal) = (*prim[iVal])[NodeIndex(i - 1, j)] /*- 0.5 * (*limiterCoef[iVal])[NodeIndex(i + 1, j, k)] * grad.dot(r)*/;
				}
				riemannSolver_->Solver(riemann_para);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
					(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
			}
		}
	}
	void Solver_NS_2D_Struct::InviscidFluxMUSCL()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_Primitive;
		auto& cons = m_Conservative;
		auto& limiterCoef = m_LimiterCoef;
		auto& res = m_Residual;
		auto& coordTrans = m_CoordTrans;
		// 起始点和终止点的编号,s: start, e: end
		int is, ie, js, je;
		grid->GetRange(is, ie, js, je);
		int iNode;
		// grid->GetNodeIndex(i, j, k)的lamda表达式
		auto NodeIndex = [&](int i, int j) {return grid->GetNodeIndex(i, j); };
		auto Prim = [&](int iVal, int i, int j) {return (*m_Primitive[iVal])[NodeIndex(i, j)]; };
		RiemannSolverPara riemann_para;
		riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
		double k = GlobalData::GetDouble("muscl_k");
		for (int j = js; j < je; j++)
		{
			for (int i = is; i < ie; i++)
			{
				iNode = NodeIndex(i, j);
				auto& jacobi = (*coordTrans[32])[iNode];
				// i direction
				riemann_para.norm(0) = (*coordTrans[16])[iNode];
				riemann_para.norm(1) = (*coordTrans[17])[iNode];
				riemann_para.norm(2) = (*coordTrans[18])[iNode];
				riemann_para.nt = (*coordTrans[19])[iNode];
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
					(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
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
					(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;

				// j direction
				riemann_para.norm(0) = (*coordTrans[20])[iNode];
				riemann_para.norm(1) = (*coordTrans[21])[iNode];
				riemann_para.norm(2) = (*coordTrans[22])[iNode];
				riemann_para.nt = (*coordTrans[23])[iNode];
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
					(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
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
					(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
			}
		}
	}
	void Solver_NS_2D_Struct::CalcLimiter()
	{
		auto& grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_Primitive;
		auto& limiterCoef = m_LimiterCoef;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
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
					maxVal = minVal = (*prim[iVal])[iNode];
					for (size_t iNeib = 0; iNeib < neighborNodeIndex.size(); ++iNeib)
					{
						maxVal = Max(maxVal, (*prim[iVal])[neighborNodeIndex[iNeib]]);
						minVal = Min(minVal, (*prim[iVal])[neighborNodeIndex[iNeib]]);
					}
					double gradx = (*primGradX[iVal])[iNode];
					double grady = (*primGradY[iVal])[iNode];
					double deltaMax = maxVal - (*prim[iVal])[iNode];
					double deltaMin = minVal - (*prim[iVal])[iNode];
					double tempCoef = LARGE_NUMBER;
					(*limiterCoef[iVal])[iNode] = LARGE_NUMBER;
					for (int iNeighbor = 0; iNeighbor < neighborNodeIndex.size(); ++iNeighbor)
					{
						auto current2Neighbor = nodeCoord[neighborNodeIndex[iNeighbor]] - nodeCoord[iNode];
						double delta2 = 0.5 * (current2Neighbor(0) * gradx + current2Neighbor(1) * grady);
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
	void Solver_NS_2D_Struct::InletBC(Boundary& bound)
	{
		FlowSolverParaPtr para = GetPara();
		int  ghostIndex = bound.GetGhostIndex();
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
		auto& cons0 = *m_Conservative[0];
		auto& cons1 = *m_Conservative[1];
		auto& cons2 = *m_Conservative[2];
		auto& cons3 = *m_Conservative[3];
		auto& cons4 = *m_Conservative[4];
		double prim_init[5];
		prim_init[0] = para->GetInflowDensity();
		prim_init[1] = para->GetInflowVelocityX();
		prim_init[2] = para->GetInflowVelocityY();
		prim_init[3] = para->GetInflowVelocityZ();
		prim_init[4] = para->GetInflowPressure();
		rho[ghostIndex] = prim_init[0];
		u[ghostIndex] = prim_init[1];
		v[ghostIndex] = prim_init[2];
		w[ghostIndex] = prim_init[3];
		p[ghostIndex] = prim_init[4];
		Prim2Cons(rho[ghostIndex], u[ghostIndex], v[ghostIndex], w[ghostIndex], p[ghostIndex],
			cons0[ghostIndex], cons1[ghostIndex], cons2[ghostIndex], cons3[ghostIndex], cons4[ghostIndex]);
	}
	void Solver_NS_2D_Struct::OutletBC(Boundary& bound)
	{
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
		auto& cons0 = *m_Conservative[0];
		auto& cons1 = *m_Conservative[1];
		auto& cons2 = *m_Conservative[2];
		auto& cons3 = *m_Conservative[3];
		auto& cons4 = *m_Conservative[4];
		int boundIndex = bound.GetIndex();
		int ghostIndex = bound.GetGhostIndex();
		rho[ghostIndex] = rho[boundIndex];
		u[ghostIndex] = u[boundIndex];
		v[ghostIndex] = v[boundIndex];
		w[ghostIndex] = w[boundIndex];
		p[ghostIndex] = p[boundIndex];
		Prim2Cons(rho[ghostIndex], u[ghostIndex], v[ghostIndex], w[ghostIndex], p[ghostIndex],
			cons0[ghostIndex], cons1[ghostIndex], cons2[ghostIndex], cons3[ghostIndex], cons4[ghostIndex]);
	}
}