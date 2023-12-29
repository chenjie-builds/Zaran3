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
			if (coordTrans.J() > 1e15 || coordTrans.J() < 0|| isnan(coordTrans.J()) || isinf(coordTrans.J()))
			{
				ZaranLog::warn("jacobi is too large or negative: {}", coordTrans.J());
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
		auto& prim = m_Primitive;
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
#pragma omp parallel for private(A, b, grad, omega, deltaVal, deltaX, deltaY, deltaZ)
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
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
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
		auto& primGradX = m_PrimGradX;
		auto& primGradY = m_PrimGradY;
		auto& primGradZ = m_PrimGradZ;
		auto& limiterCoef = m_LimiterCoef;
		auto& res = m_Residual;
		auto& coordTrans = m_CoordTrans;
		RiemannSolverPara riemann_para;
		riemann_para.gammaL = riemann_para.gammaR = 1.4;
#pragma omp parallel for private(riemann_para)
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner)
				continue;
			auto& jacobi = (*coordTrans[32])[iNode];
			// i direction
			riemann_para.norm(0) = (*coordTrans[16])[iNode];
			riemann_para.norm(1) = (*coordTrans[17])[iNode];
			riemann_para.norm(2) = (*coordTrans[18])[iNode];
			riemann_para.nt = (*coordTrans[19])[iNode];
			MidPointReconstruct(templateI[iNode][1], templateI[iNode][2], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
			MidPointReconstruct(templateI[iNode][0], templateI[iNode][1], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
			// j direction
			riemann_para.norm(0) = (*coordTrans[20])[iNode];
			riemann_para.norm(1) = (*coordTrans[21])[iNode];
			riemann_para.norm(2) = (*coordTrans[22])[iNode];
			riemann_para.nt = (*coordTrans[23])[iNode];
			MidPointReconstruct(templateJ[iNode][1], templateJ[iNode][2], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
			MidPointReconstruct(templateJ[iNode][0], templateJ[iNode][1], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
			// k direction
			riemann_para.norm(0) = (*coordTrans[24])[iNode];
			riemann_para.norm(1) = (*coordTrans[25])[iNode];
			riemann_para.norm(2) = (*coordTrans[26])[iNode];
			riemann_para.nt = (*coordTrans[27])[iNode];
			MidPointReconstruct(templateK[iNode][1], templateK[iNode][2], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
			MidPointReconstruct(templateK[iNode][0], templateK[iNode][1], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
		}
	}

	void Solver_NS_3D::ViscousFlux()
	{
	}

	void Solver_NS_3D::SourceFlux()
	{
	}



}