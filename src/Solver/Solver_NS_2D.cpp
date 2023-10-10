#include"Solver_NS_2D.h"
namespace zaran
{
	void Solver_NS_2D::InitField()
	{
		GridPtr grid = GetGrid();
		auto& rho = *m_Primtive[0];
		auto& u = *m_Primtive[1];
		auto& v = *m_Primtive[2];
		auto& w = *m_Primtive[3];
		auto& p = *m_Primtive[4];
		auto& NodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = NodeTopo->GetCoordinate();
		FlowSolverParaPtr para = GetPara();
		int initType = para->GetInitFieldType();
		DVector primInit = para->GetPrimitiveInflow();
		primInit[0] = 6.4;
		primInit[1] = 3.125;
		primInit[2] = 0;
		primInit[3] = 0;
		primInit[4] = 18.5;
		int nTotalNodeNum = grid->GetTotalNodeNum();
		double x, y, z;
		para->SetPrimitiveInflow(primInit);
		for (int iNode = 0; iNode < nTotalNodeNum; ++iNode)
		{
			x = nodeCoord[iNode].x();
			y = nodeCoord[iNode].y();
			if (x <= 0.1)
			{
				rho[iNode] = 6.45;
				u[iNode] = 3.125;
				v[iNode] = 0.0;
				w[iNode] = 0.0;
				p[iNode] = 18.5;
			}
			else
			{
				rho[iNode] = 1.4;
				u[iNode] = 0;
				v[iNode] = 0;
				w[iNode] = 0.0;
				p[iNode] = 1.0;
			}
		}
		Primitive2Conservative();
	}
	void Solver_NS_2D::ComputeCoordTrans()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeType = nodeTopo->GetType();
		auto& tempI = nodeTopo->GetTemplateI();
		auto& tempJ = nodeTopo->GetTemplateJ();
		auto& coordTransCoef = m_CoordTrans;
		int nInnerNode = grid->GetInnerNodeNum();
		CoordTrans coordTrans;
		DVector3D xRight, xLeft, yRight, yLeft;
		for (size_t iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner && nodeType[iNode] != NodeType::hole)
				continue;
			xLeft = nodeCoord[tempI[iNode][0]];
			xRight = nodeCoord[tempI[iNode][2]];
			yLeft = nodeCoord[tempJ[iNode][0]];
			yRight = nodeCoord[tempJ[iNode][2]];
			coordTrans.CalcCoordTrans(int(grid->GetDimension()), xRight, xLeft, yRight, yLeft);
			if (coordTrans.J() < 0)
			{
				tempJ[iNode] = IArray{ tempJ[iNode][2], tempJ[iNode][1], tempJ[iNode][0] };
				coordTrans.CalcCoordTrans(int(grid->GetDimension()), xRight, xLeft, yLeft, yRight);
			}
			if (abs(coordTrans.J()) < SMALL_NUMBER || isnan(abs(coordTrans.J())) || isinf((coordTrans.J())))
			{
				IArray currentTempI = { tempI[iNode][0], tempI[iNode][1], tempJ[iNode][0] };
				IArray currentTempJ = { tempI[iNode][2], tempJ[iNode][1], tempJ[iNode][0] };
				tempI[iNode] = currentTempI;
				tempJ[iNode] = currentTempJ;
				coordTrans.CalcCoordTrans(int(grid->GetDimension()), xLeft, yLeft, xRight, yRight);
				if (coordTrans.J() < 0)
				{
					tempJ[iNode] = IArray{ tempJ[iNode][2], tempJ[iNode][1], tempJ[iNode][0] };
					coordTrans.CalcCoordTrans(int(grid->GetDimension()), xLeft, yLeft, yRight, xRight);
				}
			}

			// check coordinate
			if (coordTrans.J() > 1e15 || coordTrans.J() < 0)
			{
				ZaranLog::warn("Node {}: {},{},{}", iNode, nodeCoord[iNode].x(), nodeCoord[iNode].y(), nodeCoord[iNode].z());
				ZaranLog::info("xLeft index={}: {},{},{}", tempI[iNode][0], xLeft.x(), xLeft.y(), xLeft.z());
				ZaranLog::info("xRight index={}: {},{},{}", tempI[iNode][2], xRight.x(), xRight.y(), xRight.z());
				ZaranLog::info("yLeft index={}: {},{},{}", tempJ[iNode][0], yLeft.x(), yLeft.y(), yLeft.z());
				ZaranLog::info("yRight index={}: {},{},{}", tempJ[iNode][2], yRight.x(), yRight.y(), yRight.z());
			}

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
	void Solver_NS_2D::ComputeGradientWLS()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		auto& prim = m_Primtive;
		auto& limiterCoef = m_LimiterCoef;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		int nInnerNode = grid->GetInnerNodeNum();
		int nBoundNode = grid->GetBoundNodeNum();
		Matrix2d A;
		DVector2D b, grad;
		double omega = 0;
		double deltaVal;
		double deltaX, deltaY;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			auto& currentCoord = nodeCoord[iNode];
			auto& neighborNodeVec = nodeNeighbor[iNode];
			for (size_t iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				A.setZero();
				b.setZero();
				for (size_t iNeib = 0; iNeib < neighborNodeVec.size(); ++iNeib)
				{
					omega = DistanceOfTwoPoints(nodeCoord[neighborNodeVec[iNeib]].data(), nodeCoord[iNode].data());
					if (abs(omega) < SMALL_NUMBER)
						continue;
					omega = 1.0 / omega;
					deltaVal = (*prim[iVal])[neighborNodeVec[iNeib]] - (*prim[iVal])[iNode];
					deltaX = nodeCoord[neighborNodeVec[iNeib]].x() - nodeCoord[iNode].x();
					deltaY = nodeCoord[neighborNodeVec[iNeib]].y() - nodeCoord[iNode].y();
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
	void Solver_NS_2D::ComputeTimeStepLocal()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		FlowSolverParaPtr para = GetPara();
		auto& rho = *m_Primtive[0];
		auto& u = *m_Primtive[1];
		auto& v = *m_Primtive[2];
		auto& w = *m_Primtive[3];
		auto& p = *m_Primtive[4];
		auto& dt = *m_TimeStep;
		auto& coordTrans = m_CoordTrans;
		double cfl = para->GetCflNumber();
		int nInnerNode = grid->GetInnerNodeNum();
		double minDt = LARGE_NUMBER;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{

			double gamma = 1.4;
			if (nodeType[iNode] != NodeType::inner)
				continue;
			double c = sqrt(gamma * p[iNode] / rho[iNode]);
			double normXi = sqrt((*coordTrans[16])[iNode] * (*coordTrans[16])[iNode] + (*coordTrans[17])[iNode] * (*coordTrans[17])[iNode] + (*coordTrans[18])[iNode] * (*coordTrans[19])[iNode]);
			double normEta = sqrt((*coordTrans[20])[iNode] * (*coordTrans[20])[iNode] + (*coordTrans[21])[iNode] * (*coordTrans[21])[iNode] + (*coordTrans[22])[iNode] * (*coordTrans[22])[iNode]);
			double uXi = u[iNode] * (*coordTrans[16])[iNode] + v[iNode] * (*coordTrans[17])[iNode] + w[iNode] * (*coordTrans[18])[iNode] + (*coordTrans[19])[iNode];
			double uEta = u[iNode] * (*coordTrans[20])[iNode] + v[iNode] * (*coordTrans[21])[iNode] + w[iNode] * (*coordTrans[22])[iNode] + (*coordTrans[23])[iNode];
			double lamda = abs(uXi) + abs(uEta) + c * (normXi + normEta);
			dt[iNode] = cfl / lamda;
			minDt = Min(minDt, dt[iNode]);
		}
		GlobalData::Update("dt", minDt);
	}

	void Solver_NS_2D::InviscidFlux()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& templateI = nodeTopo->GetTemplateI();
		auto& templateJ = nodeTopo->GetTemplateJ();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_Primtive;
		auto& cons = m_Conservative;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& limiterCoef = m_LimiterCoef;
		auto& res = m_Residual;
		auto& coordTrans = m_CoordTrans;
		int nInnerNode = grid->GetInnerNodeNum();
		int nBoundNode = grid->GetBoundNodeNum();
		DVector2D r, grad;
		Ptr<RiemannSolverPara >riemanPara = std::make_shared<RiemannSolverPara>();
		riemanPara->gammaL = riemanPara->gammaR = 1.4;
		for (size_t iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner)
				continue;
			auto& jacobi = (*coordTrans[32])[iNode];
			// i direction
			riemanPara->norm(0) = (*coordTrans[16])[iNode];
			riemanPara->norm(1) = (*coordTrans[17])[iNode];
			riemanPara->norm(2) = 0;
			riemanPara->nt = (*coordTrans[19])[iNode];
			auto& tempI = templateI[iNode];
			r[0] = nodeCoord[tempI[2]][0] - nodeCoord[iNode][0];
			r[1] = nodeCoord[tempI[2]][1] - nodeCoord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				riemanPara->primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempI[2]];
				grad(1) = (*primGradY[iVal])[tempI[2]];
				riemanPara->primR(iVal) = (*prim[iVal])[tempI[2]] - 0.5 * (*limiterCoef[iVal])[tempI[2]] * grad.dot(r);
			}

			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
			r[0] = nodeCoord[tempI[0]][0] - nodeCoord[iNode][0];
			r[1] = nodeCoord[tempI[0]][1] - nodeCoord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempI[0]];
				grad(1) = (*primGradY[iVal])[tempI[0]];
				riemanPara->primL(iVal) = (*prim[iVal])[tempI[0]] - 0.5 * (*limiterCoef[iVal])[tempI[0]] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				riemanPara->primR(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}

			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;

			// j direction
			riemanPara->norm(0) = (*coordTrans[20])[iNode];
			riemanPara->norm(1) = (*coordTrans[21])[iNode];
			riemanPara->norm(2) = 0;
			riemanPara->nt = (*coordTrans[23])[iNode];
			auto& tempJ = templateJ[iNode];
			r[0] = nodeCoord[tempJ[2]][0] - nodeCoord[iNode][0];
			r[1] = nodeCoord[tempJ[2]][1] - nodeCoord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				riemanPara->primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempJ[2]];
				grad(1) = (*primGradY[iVal])[tempJ[2]];
				riemanPara->primR(iVal) = (*prim[iVal])[tempJ[2]] - 0.5 * (*limiterCoef[iVal])[tempJ[2]] * grad.dot(r);
			}
			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
			r[0] = nodeCoord[tempJ[0]][0] - nodeCoord[iNode][0];
			r[1] = nodeCoord[tempJ[0]][1] - nodeCoord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempJ[0]];
				grad(1) = (*primGradY[iVal])[tempJ[0]];
				riemanPara->primL(iVal) = (*prim[iVal])[tempJ[0]] - 0.5 * (*limiterCoef[iVal])[tempJ[0]] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				riemanPara->primR(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}

			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;
				if (isnan((*res[iVal])[iNode]) || isinf((*res[iVal])[iNode]))
					ZaranLog::error("inode={},NAN in Residual!", iNode);
			}
		}
		HoleInviscidFlux();
	}
	void Solver_NS_2D::HoleInviscidFlux()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& neibor = nodeTopo->GetNeighborCloud();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_Primtive;
		auto& cons = m_Conservative;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& limiterCoef = m_LimiterCoef;
		auto& res = m_Residual;
		auto& coordTrans = m_CoordTrans;
		int nInnerNode = grid->GetInnerNodeNum();
		int nBoundNode = grid->GetBoundNodeNum();
		DVector2D r, grad;
		Ptr<RiemannSolverPara >riemanPara = std::make_shared<RiemannSolverPara>();
		riemanPara->gammaL = riemanPara->gammaR = 1.4;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (nodeType[iNode] != NodeType::hole)
				continue;
			auto& jacobi = (*coordTrans[32])[iNode];
			// i direction
			riemanPara->norm(0) = (*coordTrans[16])[iNode];
			riemanPara->norm(1) = (*coordTrans[17])[iNode];
			riemanPara->norm(2) = 0;
			riemanPara->nt = (*coordTrans[19])[iNode];
			IArray tempI = { neibor[iNode][2],iNode,neibor[iNode][0] };
			r[0] = nodeCoord[tempI[2]][0] - nodeCoord[iNode][0];
			r[1] = nodeCoord[tempI[2]][1] - nodeCoord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				riemanPara->primL(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempI[2]];
				grad(1) = (*primGradY[iVal])[tempI[2]];
				riemanPara->primR(iVal) = (*prim[iVal])[tempI[2]] /*- 0.5 * (*limiterCoef[iVal])[tempI[2]] * grad.dot(r)*/;
			}

			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
			r[0] = nodeCoord[tempI[0]][0] - nodeCoord[iNode][0];
			r[1] = nodeCoord[tempI[0]][1] - nodeCoord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempI[0]];
				grad(1) = (*primGradY[iVal])[tempI[0]];
				riemanPara->primL(iVal) = (*prim[iVal])[tempI[0]]/* - 0.5 * (*limiterCoef[iVal])[tempI[0]] * grad.dot(r)*/;
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				riemanPara->primR(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
			}

			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;

			// j direction
			riemanPara->norm(0) = (*coordTrans[20])[iNode];
			riemanPara->norm(1) = (*coordTrans[21])[iNode];
			riemanPara->norm(2) = 0;
			riemanPara->nt = (*coordTrans[23])[iNode];
			IArray tempJ = { neibor[iNode][3],iNode,neibor[iNode][1] };
			r[0] = nodeCoord[tempJ[2]][0] - nodeCoord[iNode][0];
			r[1] = nodeCoord[tempJ[2]][1] - nodeCoord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				riemanPara->primL(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempJ[2]];
				grad(1) = (*primGradY[iVal])[tempJ[2]];
				riemanPara->primR(iVal) = (*prim[iVal])[tempJ[2]]/* - 0.5 * (*limiterCoef[iVal])[tempJ[2]] * grad.dot(r)*/;
			}
			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
			r[0] = nodeCoord[tempJ[0]][0] - nodeCoord[iNode][0];
			r[1] = nodeCoord[tempJ[0]][1] - nodeCoord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempJ[0]];
				grad(1) = (*primGradY[iVal])[tempJ[0]];
				riemanPara->primL(iVal) = (*prim[iVal])[tempJ[0]] /*- 0.5 * (*limiterCoef[iVal])[tempJ[0]] * grad.dot(r)*/;
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				riemanPara->primR(iVal) = (*prim[iVal])[iNode] /*+ 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r)*/;
			}

			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;
				if (isnan((*res[iVal])[iNode]) || isinf((*res[iVal])[iNode]))
					ZaranLog::error("inode={},NAN in Residual!", iNode);
			}
		}
	}
	void Solver_NS_2D::ViscousFlux()
	{

	}

	void Solver_NS_2D::SourceFlux()
	{
	}

}