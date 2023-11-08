#include "Solver_NS_3D.h"
namespace zaran
{
	void Solver_NS_3D::ComputeCoordTrans()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeType = nodeTopo->GetType();
		auto& tempI = nodeTopo->GetTemplateI();
		auto& tempJ = nodeTopo->GetTemplateJ();
		auto& tempK = nodeTopo->GetTemplateK();
		auto& coordTransCoef = m_CoordTrans;
		int nInnerNode = grid->GetInnerNodeNum();
		CoordTrans coordTrans;
		DVector3D xRight, xLeft, yRight, yLeft, zRight, zLeft;
		for (size_t iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner)
				continue;
			xLeft = nodeCoord[tempI[iNode][0]];
			xRight = nodeCoord[tempI[iNode][2]];
			yLeft = nodeCoord[tempJ[iNode][0]];
			yRight = nodeCoord[tempJ[iNode][2]];
			zLeft = nodeCoord[tempK[iNode][0]];
			zRight = nodeCoord[tempK[iNode][2]];
			coordTrans.CalcCoordTrans(int(grid->GetDimension()), xRight, xLeft, yRight, yLeft, zRight, zLeft);
			if (coordTrans.J() < 0)
			{
				tempK[iNode] = IArray{ tempK[iNode][2], tempK[iNode][1], tempK[iNode][0] };
				coordTrans.CalcCoordTrans(3, xRight, xLeft, yRight, yLeft, zLeft, zRight);
			}
			if (abs(coordTrans.J()) < SMALL_NUMBER || isnan(abs(coordTrans.J())) || isinf((coordTrans.J())))
			{
				IArray currentTempJ = { tempJ[iNode][0], tempJ[iNode][1], tempK[iNode][0] };
				IArray currentTempK = { tempJ[iNode][2], tempK[iNode][1], tempK[iNode][2] };
				tempJ[iNode] = currentTempJ;
				tempK[iNode] = currentTempK;
				coordTrans.CalcCoordTrans(3, xRight, xLeft, yRight, zRight, yLeft, zLeft);
				if (coordTrans.J() < 0)
				{
					tempK[iNode] = IArray{ tempK[iNode][2], tempK[iNode][1], tempK[iNode][0] };
					coordTrans.CalcCoordTrans(3, xRight, xLeft, yRight, yLeft, zLeft, zRight);
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
				ZaranLog::info("zLeft index={}: {},{},{}", tempK[iNode][0], zLeft.x(), zLeft.y(), zLeft.z());
				ZaranLog::info("zRight index={}: {},{},{}", tempK[iNode][2], zRight.x(), zRight.y(), zRight.z());

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

	void Solver_NS_3D::ComputeGradientWLS()
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
		Matrix3d A;
		DVector3D b, grad;
		double omega = 0;
		double deltaVal;
		double deltaX, deltaY, deltaZ;
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
					deltaZ = nodeCoord[neighborNodeVec[iNeib]].z() - nodeCoord[iNode].z();
					A(0, 0) += omega * deltaX * deltaX;
					A(0, 1) += omega * deltaX * deltaY;
					A(0, 2) += omega * deltaX * deltaZ;
					A(1, 0) += omega * deltaY * deltaX;
					A(1, 1) += omega * deltaY * deltaY;
					A(1, 2) += omega * deltaY * deltaZ;
					A(2, 0) += omega * deltaZ * deltaX;
					A(2, 1) += omega * deltaZ * deltaY;
					A(2, 2) += omega * deltaZ * deltaZ;
					b(0) += omega * deltaVal * deltaX;
					b(1) += omega * deltaVal * deltaY;
					b(2) += omega * deltaVal * deltaZ;
				}
				grad = A.inverse() * b;
				(*primGradX[iVal])[iNode] = grad.x();
				(*primGradY[iVal])[iNode] = grad.y();
				(*primGradZ[iVal])[iNode] = grad.z();
			}
		}
	}

	void Solver_NS_3D::ComputeTimeStepLocal()
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
			double normZeta = sqrt((*coordTrans[24])[iNode] * (*coordTrans[24])[iNode] + (*coordTrans[25])[iNode] * (*coordTrans[25])[iNode] + (*coordTrans[26])[iNode] * (*coordTrans[26])[iNode]);
			double uXi = u[iNode] * (*coordTrans[16])[iNode] + v[iNode] * (*coordTrans[17])[iNode] + w[iNode] * (*coordTrans[18])[iNode] + (*coordTrans[19])[iNode];
			double uEta = u[iNode] * (*coordTrans[20])[iNode] + v[iNode] * (*coordTrans[21])[iNode] + w[iNode] * (*coordTrans[22])[iNode] + (*coordTrans[23])[iNode];
			double uZeta = u[iNode] * (*coordTrans[24])[iNode] + v[iNode] * (*coordTrans[25])[iNode] + w[iNode] * (*coordTrans[26])[iNode] + (*coordTrans[27])[iNode];
			double lamda = abs(uXi) + abs(uEta) + abs(uZeta) + c * (normXi + normEta + normZeta);
			dt[iNode] = cfl / lamda;
			minDt = Min(minDt, dt[iNode]);
		}
		GlobalData::Update("dt", minDt);
	}

	void Solver_NS_3D::InviscidFlux()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& templateI = nodeTopo->GetTemplateI();
		auto& templateJ = nodeTopo->GetTemplateJ();
		auto& templateK = nodeTopo->GetTemplateK();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_Primtive;
		auto& cons = m_Conservative;
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		auto& limiterCoef = m_LimiterCoef;
		auto& res = m_Residual;
		auto& coordTrans = m_CoordTrans;
		int nInnerNode = grid->GetInnerNodeNum();
		int nBoundNode = grid->GetBoundNodeNum();
		DVector3D r, grad;
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
			riemanPara->norm(2) = (*coordTrans[18])[iNode];
			riemanPara->nt = (*coordTrans[19])[iNode];
			auto& tempI = templateI[iNode];
			r = nodeCoord[tempI[2]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemanPara->primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempI[2]];
				grad(1) = (*primGradY[iVal])[tempI[2]];
				grad(2) = (*primGradZ[iVal])[tempI[2]];
				riemanPara->primR(iVal) = (*prim[iVal])[tempI[2]] - 0.5 * (*limiterCoef[iVal])[tempI[2]] * grad.dot(r);
			}

			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
			r = nodeCoord[tempI[0]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempI[0]];
				grad(1) = (*primGradY[iVal])[tempI[0]];
				grad(2) = (*primGradZ[iVal])[tempI[0]];
				riemanPara->primL(iVal) = (*prim[iVal])[tempI[0]] - 0.5 * (*limiterCoef[iVal])[tempI[0]] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemanPara->primR(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}

			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;

			// j direction
			riemanPara->norm(0) = (*coordTrans[20])[iNode];
			riemanPara->norm(1) = (*coordTrans[21])[iNode];
			riemanPara->norm(2) = (*coordTrans[22])[iNode];
			riemanPara->nt = (*coordTrans[23])[iNode];
			auto& tempJ = templateJ[iNode];
			r = nodeCoord[tempJ[2]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemanPara->primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempJ[2]];
				grad(1) = (*primGradY[iVal])[tempJ[2]];
				grad(2) = (*primGradZ[iVal])[tempJ[2]];
				riemanPara->primR(iVal) = (*prim[iVal])[tempJ[2]] - 0.5 * (*limiterCoef[iVal])[tempJ[2]] * grad.dot(r);
			}
			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
			r = nodeCoord[tempJ[0]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempJ[0]];
				grad(1) = (*primGradY[iVal])[tempJ[0]];
				grad(2) = (*primGradZ[iVal])[tempJ[0]];
				riemanPara->primL(iVal) = (*prim[iVal])[tempJ[0]] - 0.5 * (*limiterCoef[iVal])[tempJ[0]] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemanPara->primR(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}

			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemanPara->flux[iVal] / jacobi;

			// k direction
			riemanPara->norm(0) = (*coordTrans[24])[iNode];
			riemanPara->norm(1) = (*coordTrans[25])[iNode];
			riemanPara->norm(2) = (*coordTrans[26])[iNode];
			riemanPara->nt = (*coordTrans[27])[iNode];
			auto& tempK = templateK[iNode];
			r = nodeCoord[tempK[2]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemanPara->primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempK[2]];
				grad(1) = (*primGradY[iVal])[tempK[2]];
				grad(2) = (*primGradZ[iVal])[tempK[2]];
				riemanPara->primR(iVal) = (*prim[iVal])[tempK[2]] - 0.5 * (*limiterCoef[iVal])[tempK[2]] * grad.dot(r);
			}
			riemannSolver_->Solver(riemanPara);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemanPara->flux[iVal] / jacobi;
			r = nodeCoord[tempK[0]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempK[0]];
				grad(1) = (*primGradY[iVal])[tempK[0]];
				grad(2) = (*primGradZ[iVal])[tempK[0]];
				riemanPara->primL(iVal) = (*prim[iVal])[tempK[0]] - 0.5 * (*limiterCoef[iVal])[tempK[0]] * grad.dot(r);
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
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
	}

	void Solver_NS_3D::ViscousFlux()
	{
	}

	void Solver_NS_3D::SourceFlux()
	{
	}

}