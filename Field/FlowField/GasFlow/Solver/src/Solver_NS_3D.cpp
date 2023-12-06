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
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& prim = m_Primitive;
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
		RiemannSolverPara riemann_para;
		riemann_para.gammaL = riemann_para.gammaR = 1.4;
#pragma omp parallel for private(r, grad, riemann_para)
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
			auto& tempI = templateI[iNode];
			r = nodeCoord[tempI[2]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemann_para.primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			if (riemann_para.primL(0) <= 0 || riemann_para.primL(4) <= 0)
			{
				ZaranLog::info("neg (i+1/2)L value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primL(iVal) = (*prim[iVal])[iNode];
				}
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempI[2]];
				grad(1) = (*primGradY[iVal])[tempI[2]];
				grad(2) = (*primGradZ[iVal])[tempI[2]];
				riemann_para.primR(iVal) = (*prim[iVal])[tempI[2]] - 0.5 * (*limiterCoef[iVal])[tempI[2]] * grad.dot(r);
			}
			if (riemann_para.primR(0) <= 0 || riemann_para.primR(4) <= 0)
			{
				ZaranLog::info("neg (i+1/2)R value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primR(iVal) = (*prim[iVal])[tempI[2]];
				}
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;

			r = nodeCoord[tempI[0]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempI[0]];
				grad(1) = (*primGradY[iVal])[tempI[0]];
				grad(2) = (*primGradZ[iVal])[tempI[0]];
				riemann_para.primL(iVal) = (*prim[iVal])[tempI[0]] - 0.5 * (*limiterCoef[iVal])[tempI[0]] * grad.dot(r);
			}
			if (riemann_para.primL(0) <= 0 || riemann_para.primL(4) <= 0)
			{
				ZaranLog::info("neg (i-1/2)L value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primL(iVal) = (*prim[iVal])[tempI[0]];
				}
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemann_para.primR(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			if (riemann_para.primR(0) <= 0 || riemann_para.primR(4) <= 0)
			{
				ZaranLog::info("neg (i-1/2)R value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primR(iVal) = (*prim[iVal])[iNode];
				}
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
			// j direction
			riemann_para.norm(0) = (*coordTrans[20])[iNode];
			riemann_para.norm(1) = (*coordTrans[21])[iNode];
			riemann_para.norm(2) = (*coordTrans[22])[iNode];
			riemann_para.nt = (*coordTrans[23])[iNode];
			auto& tempJ = templateJ[iNode];
			r = nodeCoord[tempJ[2]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemann_para.primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			if (riemann_para.primL(0) <= 0 || riemann_para.primL(4) <= 0)
			{
				ZaranLog::info("neg (j+1/2)L value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primL(iVal) = (*prim[iVal])[iNode];
				}
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempJ[2]];
				grad(1) = (*primGradY[iVal])[tempJ[2]];
				grad(2) = (*primGradZ[iVal])[tempJ[2]];
				riemann_para.primR(iVal) = (*prim[iVal])[tempJ[2]] - 0.5 * (*limiterCoef[iVal])[tempJ[2]] * grad.dot(r);
			}
			if (riemann_para.primR(0) <= 0 || riemann_para.primR(4) <= 0)
			{
				ZaranLog::info("neg (j+1/2)R value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primR(iVal) = (*prim[iVal])[tempJ[2]];
				}
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
			r = nodeCoord[tempJ[0]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempJ[0]];
				grad(1) = (*primGradY[iVal])[tempJ[0]];
				grad(2) = (*primGradZ[iVal])[tempJ[0]];
				riemann_para.primL(iVal) = (*prim[iVal])[tempJ[0]] - 0.5 * (*limiterCoef[iVal])[tempJ[0]] * grad.dot(r);
			}
			if (riemann_para.primL(0) <= 0 || riemann_para.primL(4) <= 0)
			{
				ZaranLog::info("neg (j-1/2)L value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primL(iVal) = (*prim[iVal])[tempJ[0]];
				}
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemann_para.primR(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			if (riemann_para.primR(0) <= 0 || riemann_para.primR(4) <= 0)
			{
				ZaranLog::info("neg (j-1/2)R value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primR(iVal) = (*prim[iVal])[iNode];
				}
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
			// k direction
			riemann_para.norm(0) = (*coordTrans[24])[iNode];
			riemann_para.norm(1) = (*coordTrans[25])[iNode];
			riemann_para.norm(2) = (*coordTrans[26])[iNode];
			riemann_para.nt = (*coordTrans[27])[iNode];
			auto& tempK = templateK[iNode];
			r = nodeCoord[tempK[2]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemann_para.primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			if (riemann_para.primL(0) <= 0 || riemann_para.primL(4) <= 0)
			{
				ZaranLog::info("neg (k+1/2)L value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primL(iVal) = (*prim[iVal])[iNode];
				}
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempK[2]];
				grad(1) = (*primGradY[iVal])[tempK[2]];
				grad(2) = (*primGradZ[iVal])[tempK[2]];
				riemann_para.primR(iVal) = (*prim[iVal])[tempK[2]] - 0.5 * (*limiterCoef[iVal])[tempK[2]] * grad.dot(r);
			}
			if (riemann_para.primR(0) <= 0 || riemann_para.primR(4) <= 0)
			{
				ZaranLog::info("neg (k+1/2)R value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primR(iVal) = (*prim[iVal])[tempK[2]];
				}
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
			r = nodeCoord[tempK[0]] - nodeCoord[iNode];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[tempK[0]];
				grad(1) = (*primGradY[iVal])[tempK[0]];
				grad(2) = (*primGradZ[iVal])[tempK[0]];
				riemann_para.primL(iVal) = (*prim[iVal])[tempK[0]] - 0.5 * (*limiterCoef[iVal])[tempK[0]] * grad.dot(r);
			}
			if (riemann_para.primL(0) <= 0 || riemann_para.primL(4) <= 0)
			{
				ZaranLog::info("neg (k-1/2)L value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primL(iVal) = (*prim[iVal])[tempK[0]];
				}
			}
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				grad(0) = (*primGradX[iVal])[iNode];
				grad(1) = (*primGradY[iVal])[iNode];
				grad(2) = (*primGradZ[iVal])[iNode];
				riemann_para.primR(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiterCoef[iVal])[iNode] * grad.dot(r);
			}
			if (riemann_para.primR(0) <= 0 || riemann_para.primR(4) <= 0)
			{
				ZaranLog::info("neg (k-1/2)R value in node: {}", iNode);
				for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				{
					riemann_para.primR(iVal) = (*prim[iVal])[iNode];
				}
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
				if (isnan((*res[iVal])[iNode]) || isinf((*res[iVal])[iNode]))
				{
					ZaranLog::error("step={}", GlobalData::GetInt("step"));
					ZaranLog::error("inode={},NAN in Residual!", iNode);
					ZaranLog::error("neighbor={},{},{},{},{},{}", tempI[0], tempI[2], tempJ[0], tempJ[2], tempK[0], tempK[2]);
					ZaranLog::error("res={},{},{},{},{}", (*res[0])[iNode], (*res[1])[iNode], (*res[2])[iNode], (*res[3])[iNode], (*res[4])[iNode]);
					ZaranLog::error("coord={},{},{}", nodeCoord[iNode].x(), nodeCoord[iNode].y(), nodeCoord[iNode].z());
					ZaranLog::error("coordI0={},{},{}", nodeCoord[tempI[0]].x(), nodeCoord[tempI[0]].y(), nodeCoord[tempI[0]].z());
					ZaranLog::error("coordI2={},{},{}", nodeCoord[tempI[2]].x(), nodeCoord[tempI[2]].y(), nodeCoord[tempI[2]].z());
					ZaranLog::error("coordJ0={},{},{}", nodeCoord[tempJ[0]].x(), nodeCoord[tempJ[0]].y(), nodeCoord[tempJ[0]].z());
					ZaranLog::error("coordJ2={},{},{}", nodeCoord[tempJ[2]].x(), nodeCoord[tempJ[2]].y(), nodeCoord[tempJ[2]].z());
					ZaranLog::error("coordK0={},{},{}", nodeCoord[tempK[0]].x(), nodeCoord[tempK[0]].y(), nodeCoord[tempK[0]].z());
					ZaranLog::error("coordK2={},{},{}", nodeCoord[tempK[2]].x(), nodeCoord[tempK[2]].y(), nodeCoord[tempK[2]].z());
					ZaranLog::error("node type={}", (int)nodeType[iNode]);
					ZaranLog::error("node typeI0={}", (int)nodeType[tempI[0]]);
					ZaranLog::error("node typeI2={}", (int)nodeType[tempI[2]]);
					ZaranLog::error("node typeJ0={}", (int)nodeType[tempJ[0]]);
					ZaranLog::error("node typeJ2={}", (int)nodeType[tempJ[2]]);
					ZaranLog::error("node typeK0={}", (int)nodeType[tempK[0]]);
					ZaranLog::error("node typeK2={}", (int)nodeType[tempK[2]]);
					ZaranLog::error("primI0={},{},{},{},{}", (*prim[0])[tempI[0]], (*prim[1])[tempI[0]], (*prim[2])[tempI[0]], (*prim[3])[tempI[0]], (*prim[4])[tempI[0]]);
					ZaranLog::error("primI1={},{},{},{},{}", (*prim[0])[iNode], (*prim[1])[iNode], (*prim[2])[iNode], (*prim[3])[iNode], (*prim[4])[iNode]);
					ZaranLog::error("primI2={},{},{},{},{}", (*prim[0])[tempI[2]], (*prim[1])[tempI[2]], (*prim[2])[tempI[2]], (*prim[3])[tempI[2]], (*prim[4])[tempI[2]]);
					ZaranLog::error("primJ0={},{},{},{},{}", (*prim[0])[tempJ[0]], (*prim[1])[tempJ[0]], (*prim[2])[tempJ[0]], (*prim[3])[tempJ[0]], (*prim[4])[tempJ[0]]);
					ZaranLog::error("primJ1={},{},{},{},{}", (*prim[0])[iNode], (*prim[1])[iNode], (*prim[2])[iNode], (*prim[3])[iNode], (*prim[4])[iNode]);
					ZaranLog::error("primJ2={},{},{},{},{}", (*prim[0])[tempJ[2]], (*prim[1])[tempJ[2]], (*prim[2])[tempJ[2]], (*prim[3])[tempJ[2]], (*prim[4])[tempJ[2]]);
					ZaranLog::error("primK0={},{},{},{},{}", (*prim[0])[tempK[0]], (*prim[1])[tempK[0]], (*prim[2])[tempK[0]], (*prim[3])[tempK[0]], (*prim[4])[tempK[0]]);
					ZaranLog::error("primK1={},{},{},{},{}", (*prim[0])[iNode], (*prim[1])[iNode], (*prim[2])[iNode], (*prim[3])[iNode], (*prim[4])[iNode]);
					ZaranLog::error("primK2={},{},{},{},{}", (*prim[0])[tempK[2]], (*prim[1])[tempK[2]], (*prim[2])[tempK[2]], (*prim[3])[tempK[2]], (*prim[4])[tempK[2]]);
					ZaranLog::error("limiterI0={},{},{},{},{}", (*limiterCoef[0])[tempI[0]], (*limiterCoef[1])[tempI[0]], (*limiterCoef[2])[tempI[0]], (*limiterCoef[3])[tempI[0]], (*limiterCoef[4])[tempI[0]]);
					ZaranLog::error("limiterI1={},{},{},{},{}", (*limiterCoef[0])[iNode], (*limiterCoef[1])[iNode], (*limiterCoef[2])[iNode], (*limiterCoef[3])[iNode], (*limiterCoef[4])[iNode]);
					ZaranLog::error("limiterI2={},{},{},{},{}", (*limiterCoef[0])[tempI[2]], (*limiterCoef[1])[tempI[2]], (*limiterCoef[2])[tempI[2]], (*limiterCoef[3])[tempI[2]], (*limiterCoef[4])[tempI[2]]);
					ZaranLog::error("limiterJ0={},{},{},{},{}", (*limiterCoef[0])[tempJ[0]], (*limiterCoef[1])[tempJ[0]], (*limiterCoef[2])[tempJ[0]], (*limiterCoef[3])[tempJ[0]], (*limiterCoef[4])[tempJ[0]]);
					ZaranLog::error("limiterJ1={},{},{},{},{}", (*limiterCoef[0])[iNode], (*limiterCoef[1])[iNode], (*limiterCoef[2])[iNode], (*limiterCoef[3])[iNode], (*limiterCoef[4])[iNode]);
					ZaranLog::error("limiterJ2={},{},{},{},{}", (*limiterCoef[0])[tempJ[2]], (*limiterCoef[1])[tempJ[2]], (*limiterCoef[2])[tempJ[2]], (*limiterCoef[3])[tempJ[2]], (*limiterCoef[4])[tempJ[2]]);
					ZaranLog::error("limiterK0={},{},{},{},{}", (*limiterCoef[0])[tempK[0]], (*limiterCoef[1])[tempK[0]], (*limiterCoef[2])[tempK[0]], (*limiterCoef[3])[tempK[0]], (*limiterCoef[4])[tempK[0]]);
					ZaranLog::error("limiterK1={},{},{},{},{}", (*limiterCoef[0])[iNode], (*limiterCoef[1])[iNode], (*limiterCoef[2])[iNode], (*limiterCoef[3])[iNode], (*limiterCoef[4])[iNode]);
					ZaranLog::error("limiterK2={},{},{},{},{}", (*limiterCoef[0])[tempK[2]], (*limiterCoef[1])[tempK[2]], (*limiterCoef[2])[tempK[2]], (*limiterCoef[3])[tempK[2]], (*limiterCoef[4])[tempK[2]]);
					ZaranLog::error("rhoGradI0={},{},{},{}", (*primGradX[0])[tempI[0]], (*primGradY[0])[tempI[0]], (*primGradZ[0])[tempI[0]], (*limiterCoef[0])[tempI[0]]);
					ZaranLog::error("rhoGradI1={},{},{},{}", (*primGradX[0])[iNode], (*primGradY[0])[iNode], (*primGradZ[0])[iNode], (*limiterCoef[0])[iNode]);
					ZaranLog::error("rhoGradI2={},{},{},{}", (*primGradX[0])[tempI[2]], (*primGradY[0])[tempI[2]], (*primGradZ[0])[tempI[2]], (*limiterCoef[0])[tempI[2]]);
					ZaranLog::error("rhoGradJ0={},{},{},{}", (*primGradX[0])[tempJ[0]], (*primGradY[0])[tempJ[0]], (*primGradZ[0])[tempJ[0]], (*limiterCoef[0])[tempJ[0]]);
					ZaranLog::error("rhoGradJ1={},{},{},{}", (*primGradX[0])[iNode], (*primGradY[0])[iNode], (*primGradZ[0])[iNode], (*limiterCoef[0])[iNode]);
					ZaranLog::error("rhoGradJ2={},{},{},{}", (*primGradX[0])[tempJ[2]], (*primGradY[0])[tempJ[2]], (*primGradZ[0])[tempJ[2]], (*limiterCoef[0])[tempJ[2]]);
					ZaranLog::error("rhoGradK0={},{},{},{}", (*primGradX[0])[tempK[0]], (*primGradY[0])[tempK[0]], (*primGradZ[0])[tempK[0]], (*limiterCoef[0])[tempK[0]]);
					ZaranLog::error("rhoGradK1={},{},{},{}", (*primGradX[0])[iNode], (*primGradY[0])[iNode], (*primGradZ[0])[iNode], (*limiterCoef[0])[iNode]);
					ZaranLog::error("rhoGradK2={},{},{},{}", (*primGradX[0])[tempK[2]], (*primGradY[0])[tempK[2]], (*primGradZ[0])[tempK[2]], (*limiterCoef[0])[tempK[2]]);
					ZaranLog::error("uGradI0={},{},{},{}", (*primGradX[1])[tempI[0]], (*primGradY[1])[tempI[0]], (*primGradZ[1])[tempI[0]], (*limiterCoef[1])[tempI[0]]);
					ZaranLog::error("uGradI1={},{},{},{}", (*primGradX[1])[iNode], (*primGradY[1])[iNode], (*primGradZ[1])[iNode], (*limiterCoef[1])[iNode]);
					ZaranLog::error("uGradI2={},{},{},{}", (*primGradX[1])[tempI[2]], (*primGradY[1])[tempI[2]], (*primGradZ[1])[tempI[2]], (*limiterCoef[1])[tempI[2]]);
					ZaranLog::error("uGradJ0={},{},{},{}", (*primGradX[1])[tempJ[0]], (*primGradY[1])[tempJ[0]], (*primGradZ[1])[tempJ[0]], (*limiterCoef[1])[tempJ[0]]);
					ZaranLog::error("uGradJ1={},{},{},{}", (*primGradX[1])[iNode], (*primGradY[1])[iNode], (*primGradZ[1])[iNode], (*limiterCoef[1])[iNode]);
					ZaranLog::error("uGradJ2={},{},{},{}", (*primGradX[1])[tempJ[2]], (*primGradY[1])[tempJ[2]], (*primGradZ[1])[tempJ[2]], (*limiterCoef[1])[tempJ[2]]);
					ZaranLog::error("uGradK0={},{},{},{}", (*primGradX[1])[tempK[0]], (*primGradY[1])[tempK[0]], (*primGradZ[1])[tempK[0]], (*limiterCoef[1])[tempK[0]]);
					ZaranLog::error("uGradK1={},{},{},{}", (*primGradX[1])[iNode], (*primGradY[1])[iNode], (*primGradZ[1])[iNode], (*limiterCoef[1])[iNode]);
					ZaranLog::error("uGradK2={},{},{},{}", (*primGradX[1])[tempK[2]], (*primGradY[1])[tempK[2]], (*primGradZ[1])[tempK[2]], (*limiterCoef[1])[tempK[2]]);
					ZaranLog::error("vGradI0={},{},{},{}", (*primGradX[2])[tempI[0]], (*primGradY[2])[tempI[0]], (*primGradZ[2])[tempI[0]], (*limiterCoef[2])[tempI[0]]);
					ZaranLog::error("vGradI1={},{},{},{}", (*primGradX[2])[iNode], (*primGradY[2])[iNode], (*primGradZ[2])[iNode], (*limiterCoef[2])[iNode]);
					ZaranLog::error("vGradI2={},{},{},{}", (*primGradX[2])[tempI[2]], (*primGradY[2])[tempI[2]], (*primGradZ[2])[tempI[2]], (*limiterCoef[2])[tempI[2]]);
					ZaranLog::error("vGradJ0={},{},{},{}", (*primGradX[2])[tempJ[0]], (*primGradY[2])[tempJ[0]], (*primGradZ[2])[tempJ[0]], (*limiterCoef[2])[tempJ[0]]);
					ZaranLog::error("vGradJ1={},{},{},{}", (*primGradX[2])[iNode], (*primGradY[2])[iNode], (*primGradZ[2])[iNode], (*limiterCoef[2])[iNode]);
					ZaranLog::error("vGradJ2={},{},{},{}", (*primGradX[2])[tempJ[2]], (*primGradY[2])[tempJ[2]], (*primGradZ[2])[tempJ[2]], (*limiterCoef[2])[tempJ[2]]);
					ZaranLog::error("vGradK0={},{},{},{}", (*primGradX[2])[tempK[0]], (*primGradY[2])[tempK[0]], (*primGradZ[2])[tempK[0]], (*limiterCoef[2])[tempK[0]]);
					ZaranLog::error("vGradK1={},{},{},{}", (*primGradX[2])[iNode], (*primGradY[2])[iNode], (*primGradZ[2])[iNode], (*limiterCoef[2])[iNode]);
					ZaranLog::error("vGradK2={},{},{},{}", (*primGradX[2])[tempK[2]], (*primGradY[2])[tempK[2]], (*primGradZ[2])[tempK[2]], (*limiterCoef[2])[tempK[2]]);
					ZaranLog::error("wGradI0={},{},{},{}", (*primGradX[3])[tempI[0]], (*primGradY[3])[tempI[0]], (*primGradZ[3])[tempI[0]], (*limiterCoef[3])[tempI[0]]);
					ZaranLog::error("wGradI1={},{},{},{}", (*primGradX[3])[iNode], (*primGradY[3])[iNode], (*primGradZ[3])[iNode], (*limiterCoef[3])[iNode]);
					ZaranLog::error("wGradI2={},{},{},{}", (*primGradX[3])[tempI[2]], (*primGradY[3])[tempI[2]], (*primGradZ[3])[tempI[2]], (*limiterCoef[3])[tempI[2]]);
					ZaranLog::error("wGradJ0={},{},{},{}", (*primGradX[3])[tempJ[0]], (*primGradY[3])[tempJ[0]], (*primGradZ[3])[tempJ[0]], (*limiterCoef[3])[tempJ[0]]);
					ZaranLog::error("wGradJ1={},{},{},{}", (*primGradX[3])[iNode], (*primGradY[3])[iNode], (*primGradZ[3])[iNode], (*limiterCoef[3])[iNode]);
					ZaranLog::error("wGradJ2={},{},{},{}", (*primGradX[3])[tempJ[2]], (*primGradY[3])[tempJ[2]], (*primGradZ[3])[tempJ[2]], (*limiterCoef[3])[tempJ[2]]);
					ZaranLog::error("wGradK0={},{},{},{}", (*primGradX[3])[tempK[0]], (*primGradY[3])[tempK[0]], (*primGradZ[3])[tempK[0]], (*limiterCoef[3])[tempK[0]]);
					ZaranLog::error("wGradK1={},{},{},{}", (*primGradX[3])[iNode], (*primGradY[3])[iNode], (*primGradZ[3])[iNode], (*limiterCoef[3])[iNode]);
					ZaranLog::error("wGradK2={},{},{},{}", (*primGradX[3])[tempK[2]], (*primGradY[3])[tempK[2]], (*primGradZ[3])[tempK[2]], (*limiterCoef[3])[tempK[2]]);
					ZaranLog::error("pGradI0={},{},{},{}", (*primGradX[4])[tempI[0]], (*primGradY[4])[tempI[0]], (*primGradZ[4])[tempI[0]], (*limiterCoef[4])[tempI[0]]);
					ZaranLog::error("pGradI1={},{},{},{}", (*primGradX[4])[iNode], (*primGradY[4])[iNode], (*primGradZ[4])[iNode], (*limiterCoef[4])[iNode]);
					ZaranLog::error("pGradI2={},{},{},{}", (*primGradX[4])[tempI[2]], (*primGradY[4])[tempI[2]], (*primGradZ[4])[tempI[2]], (*limiterCoef[4])[tempI[2]]);
					ZaranLog::error("pGradJ0={},{},{},{}", (*primGradX[4])[tempJ[0]], (*primGradY[4])[tempJ[0]], (*primGradZ[4])[tempJ[0]], (*limiterCoef[4])[tempJ[0]]);
					ZaranLog::error("pGradJ1={},{},{},{}", (*primGradX[4])[iNode], (*primGradY[4])[iNode], (*primGradZ[4])[iNode], (*limiterCoef[4])[iNode]);
					ZaranLog::error("pGradJ2={},{},{},{}", (*primGradX[4])[tempJ[2]], (*primGradY[4])[tempJ[2]], (*primGradZ[4])[tempJ[2]], (*limiterCoef[4])[tempJ[2]]);
					ZaranLog::error("pGradK0={},{},{},{}", (*primGradX[4])[tempK[0]], (*primGradY[4])[tempK[0]], (*primGradZ[4])[tempK[0]], (*limiterCoef[4])[tempK[0]]);
					ZaranLog::error("pGradK1={},{},{},{}", (*primGradX[4])[iNode], (*primGradY[4])[iNode], (*primGradZ[4])[iNode], (*limiterCoef[4])[iNode]);
					ZaranLog::error("pGradK2={},{},{},{}", (*primGradX[4])[tempK[2]], (*primGradY[4])[tempK[2]], (*primGradZ[4])[tempK[2]], (*limiterCoef[4])[tempK[2]]);
					exit(-1);

				}
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