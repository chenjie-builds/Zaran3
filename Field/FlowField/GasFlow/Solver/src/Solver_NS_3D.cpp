#include "Solver_NS_3D.h"
namespace zaran
{
	void Solver_NS_3D::CalcMetric()
	{
		Log::info("Compute NS 3D Coordination Transformation Coefficients");
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeType = nodeTopo->GetType();
		auto& tempI = nodeTopo->GetTemplateI();
		auto& tempJ = nodeTopo->GetTemplateJ();
		auto& tempK = nodeTopo->GetTemplateK();
		int nInnerNode = grid->GetInnerNodeNum();
		CoordTrans coordTrans;
		DVector3D xRight, xLeft, yRight, yLeft, zRight, zLeft;
		double max_jacobi, min_jacobi;
		max_jacobi = -LARGE_NUMBER;
		min_jacobi = LARGE_NUMBER;
		int max_jacobi_index = 0;
		int min_jacobi_index = 0;
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
			if (coordTrans.J() > 1e15 || coordTrans.J() < 0 || isnan(coordTrans.J()) || isinf(coordTrans.J()))
			{
				Log::warn("jacobi is too large or negative: {}", coordTrans.J());
				Log::warn("Node {}: {},{},{}", iNode, nodeCoord[iNode].x(), nodeCoord[iNode].y(), nodeCoord[iNode].z());
				Log::info("xLeft index={}: {},{},{}", tempI[iNode][0], xLeft.x(), xLeft.y(), xLeft.z());
				Log::info("xRight index={}: {},{},{}", tempI[iNode][2], xRight.x(), xRight.y(), xRight.z());
				Log::info("yLeft index={}: {},{},{}", tempJ[iNode][0], yLeft.x(), yLeft.y(), yLeft.z());
				Log::info("yRight index={}: {},{},{}", tempJ[iNode][2], yRight.x(), yRight.y(), yRight.z());
				Log::info("zLeft index={}: {},{},{}", tempK[iNode][0], zLeft.x(), zLeft.y(), zLeft.z());
				Log::info("zRight index={}: {},{},{}", tempK[iNode][2], zRight.x(), zRight.y(), zRight.z());
			}
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
			if (coordTrans.J() > max_jacobi)
			{
				max_jacobi = coordTrans.J();
				max_jacobi_index = iNode;
			}
			if (coordTrans.J() < min_jacobi)
			{
				min_jacobi = coordTrans.J();
				min_jacobi_index = iNode;
			}
		}
		Log::info("max jacobi: {}, index: {}", max_jacobi, max_jacobi_index);
		Log::info("min jacobi: {}, index: {}", min_jacobi, min_jacobi_index);
		Log::info("NS 3D Coordination Transformation Coefficients are computed");
	}

	void Solver_NS_3D::CalcGradWLS()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeCoord = nodeTopo->GetCoordinate();
		auto& nodeNeighbor = nodeTopo->GetNeighborCloud();
		int nInnerNode = grid->GetInnerNodeNum();
		int nBoundNode = grid->GetBoundNodeNum();
		Matrix3d A, A_inv;
		DVector3D b, grad;
		double omega = 0;
		double deltaVal;
		double deltaX, deltaY, deltaZ;
#pragma omp parallel for private(A,A_inv, b, grad, omega, deltaVal, deltaX, deltaY, deltaZ)
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			auto& currentCoord = nodeCoord[iNode];
			auto& neighborNodeVec = nodeNeighbor[iNode];
			A.setZero();
			for (size_t iNeib = 0; iNeib < neighborNodeVec.size(); ++iNeib)
			{
				omega = DistanceOfTwoPoints(nodeCoord[neighborNodeVec[iNeib]].data(), nodeCoord[iNode].data());
				if (abs(omega) < SMALL_NUMBER)
					continue;
				omega = 1.0 / omega;
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
			}
			A_inv = A.inverse();
			for (size_t iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				b.setZero();
				for (size_t iNeib = 0; iNeib < neighborNodeVec.size(); ++iNeib)
				{
					omega = DistanceOfTwoPoints(nodeCoord[neighborNodeVec[iNeib]].data(), nodeCoord[iNode].data());
					if (abs(omega) < SMALL_NUMBER)
						continue;
					omega = 1.0 / omega;
					deltaVal = m_prim[iVal][neighborNodeVec[iNeib]] - m_prim[iVal][iNode];
					deltaX = nodeCoord[neighborNodeVec[iNeib]].x() - nodeCoord[iNode].x();
					deltaY = nodeCoord[neighborNodeVec[iNeib]].y() - nodeCoord[iNode].y();
					deltaZ = nodeCoord[neighborNodeVec[iNeib]].z() - nodeCoord[iNode].z();
					b(0) += omega * deltaVal * deltaX;
					b(1) += omega * deltaVal * deltaY;
					b(2) += omega * deltaVal * deltaZ;
				}
				grad = A_inv * b;
				m_PrimGradX[iVal][iNode] = grad.x();
				m_PrimGradY[iVal][iNode] = grad.y();
				m_PrimGradZ[iVal][iNode] = grad.z();
			}
		}
	}

	void Solver_NS_3D::CalcTimeStepLocal()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		FlowSolverPara* para = GetPara();
		double cfl = para->GetCflNumber();
		int nInnerNode = grid->GetInnerNodeNum();
		double min_dt = LARGE_NUMBER;
#pragma omp parallel for reduction(min:min_dt)
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{

			double gamma = 1.4;
			if (nodeType[iNode] != NodeType::inner)
				continue;
			double c = sqrt(gamma * m_prim[4][iNode] / m_prim[0][iNode]);
			double normXi = sqrt(m_metric[16][iNode] * m_metric[16][iNode] + m_metric[17][iNode] * m_metric[17][iNode] + m_metric[18][iNode] * m_metric[18][iNode]);
			double normEta = sqrt(m_metric[20][iNode] * m_metric[20][iNode] + m_metric[21][iNode] * m_metric[21][iNode] + m_metric[22][iNode] * m_metric[22][iNode]);
			double normZeta = sqrt(m_metric[24][iNode] * m_metric[24][iNode] + m_metric[25][iNode] * m_metric[25][iNode] + m_metric[26][iNode] * m_metric[26][iNode]);
			double uXi = m_prim[1][iNode] * m_metric[16][iNode] + m_prim[2][iNode] * m_metric[17][iNode] + m_prim[3][iNode] * m_metric[18][iNode];
			double uEta = m_prim[1][iNode] * m_metric[20][iNode] + m_prim[2][iNode] * m_metric[21][iNode] + m_prim[3][iNode] * m_metric[22][iNode];
			double uZeta = m_prim[1][iNode] * m_metric[24][iNode] + m_prim[2][iNode] * m_metric[25][iNode] + m_prim[3][iNode] * m_metric[26][iNode];
			double lamda = abs(uXi) + abs(uEta) + abs(uZeta) + c * (normXi + normEta + normZeta);
			m_dt[iNode] = cfl / lamda;
			if (m_dt[iNode] < min_dt)
			{
				min_dt = m_dt[iNode];
			}

		}
		GlobalData::Update("dt", min_dt);
	}

	void Solver_NS_3D::InviscidFlux()
	{
		GridPtr grid = GetGrid();
		auto& nodeTopo = grid->GetNodeTopo();
		auto& nodeType = nodeTopo->GetType();
		auto& templateI = nodeTopo->GetTemplateI();
		auto& templateJ = nodeTopo->GetTemplateJ();
		auto& templateK = nodeTopo->GetTemplateK();
		RiemannSolverPara riemann_para[6];
		for (int i = 0; i < 6; ++i)
		{
			riemann_para[i].gamma_left = riemann_para[i].gamma_right = 1.4;
		}
		bool exist_negative = false;
#pragma omp parallel for private(riemann_para,exist_negative)
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (nodeType[iNode] != NodeType::inner)
				continue;
			exist_negative = false;
			auto& jacobi = m_metric[32][iNode];
			// i direction
			riemann_para[0].norm(0) = m_metric[16][iNode];
			riemann_para[0].norm(1) = m_metric[17][iNode];
			riemann_para[0].norm(2) = m_metric[18][iNode];
			riemann_para[0].nt = m_metric[19][iNode];
			riemann_para[1].norm(0) = m_metric[16][iNode];
			riemann_para[1].norm(1) = m_metric[17][iNode];
			riemann_para[1].norm(2) = m_metric[18][iNode];
			riemann_para[1].nt = m_metric[19][iNode];
			MidPointReconstruct(templateI[iNode][1], templateI[iNode][2], &riemann_para[0].prim_left(0), &riemann_para[0].prim_right(0));
			MidPointReconstruct(templateI[iNode][0], templateI[iNode][1], &riemann_para[1].prim_left(0), &riemann_para[1].prim_right(0));
			// j direction
			riemann_para[2].norm(0) = m_metric[20][iNode];
			riemann_para[2].norm(1) = m_metric[21][iNode];
			riemann_para[2].norm(2) = m_metric[22][iNode];
			riemann_para[2].nt = m_metric[23][iNode];
			riemann_para[3].norm(0) = m_metric[20][iNode];
			riemann_para[3].norm(1) = m_metric[21][iNode];
			riemann_para[3].norm(2) = m_metric[22][iNode];
			riemann_para[3].nt = m_metric[23][iNode];
			MidPointReconstruct(templateJ[iNode][1], templateJ[iNode][2], &riemann_para[2].prim_left(0), &riemann_para[2].prim_right(0));
			MidPointReconstruct(templateJ[iNode][0], templateJ[iNode][1], &riemann_para[3].prim_left(0), &riemann_para[3].prim_right(0));
			// k direction
			riemann_para[4].norm(0) = m_metric[24][iNode];
			riemann_para[4].norm(1) = m_metric[25][iNode];
			riemann_para[4].norm(2) = m_metric[26][iNode];
			riemann_para[4].nt = m_metric[27][iNode];
			riemann_para[5].norm(0) = m_metric[24][iNode];
			riemann_para[5].norm(1) = m_metric[25][iNode];
			riemann_para[5].norm(2) = m_metric[26][iNode];
			riemann_para[5].nt = m_metric[27][iNode];
			MidPointReconstruct(templateK[iNode][1], templateK[iNode][2], &riemann_para[4].prim_left(0), &riemann_para[4].prim_right(0));
			MidPointReconstruct(templateK[iNode][0], templateK[iNode][1], &riemann_para[5].prim_left(0), &riemann_para[5].prim_right(0));
			// check negative density and pressure
			for (int i = 0; i < 6; ++i)
			{
				if (riemann_para[i].prim_left(0) < 0 || riemann_para[i].prim_right(0) < 0 || riemann_para[i].prim_left(4) < 0 || riemann_para[i].prim_right(4) < 0)
				{
					exist_negative = true;
					break;
				}
			}
			if (exist_negative || m_non_physical[iNode] > 0)
			{
				MidPointReconstructFirstOrder(templateI[iNode][1], templateI[iNode][2], &riemann_para[0].prim_left(0), &riemann_para[0].prim_right(0));
				MidPointReconstructFirstOrder(templateI[iNode][0], templateI[iNode][1], &riemann_para[1].prim_left(0), &riemann_para[1].prim_right(0));
				MidPointReconstructFirstOrder(templateJ[iNode][1], templateJ[iNode][2], &riemann_para[2].prim_left(0), &riemann_para[2].prim_right(0));
				MidPointReconstructFirstOrder(templateJ[iNode][0], templateJ[iNode][1], &riemann_para[3].prim_left(0), &riemann_para[3].prim_right(0));
				MidPointReconstructFirstOrder(templateK[iNode][1], templateK[iNode][2], &riemann_para[4].prim_left(0), &riemann_para[4].prim_right(0));
				MidPointReconstructFirstOrder(templateK[iNode][0], templateK[iNode][1], &riemann_para[5].prim_left(0), &riemann_para[5].prim_right(0));
			}
			// calculate flux
			for (int i = 0; i < 6; ++i)
			{
				riemannSolver_->Solver(riemann_para[i]);
			}
			for (int iVar = 0;iVar < GetNumberOfEquations();++iVar)
			{
				m_residual[iVar][iNode] = (riemann_para[0].flux[iVar] - riemann_para[1].flux[iVar] + riemann_para[2].flux[iVar] - riemann_para[3].flux[iVar] + riemann_para[4].flux[iVar] - riemann_para[5].flux[iVar]) / jacobi;
			}
		}
	}

	void Solver_NS_3D::ViscousFlux()
	{
	}

	void Solver_NS_3D::SourceFlux()
	{
	}

	void Solver_NS_3D::CalcForce()
	{
		return;
		GridPtr grid = GetGrid();
		auto& faceTopo = grid->GetFaceTopo();
		int nFace = faceTopo->GetFaceNum();
		double face_pressure;
		double force[3] = { 0 };
		for (int iFace = 0;iFace < nFace;++iFace)
		{
			face_pressure = 0;
			int* face2node = faceTopo->GetFace2Node(iFace);
			for (int iNode = 0;iNode < faceTopo->GetFaceNodeNum(iFace);++iNode)
			{
				face_pressure += m_prim[4][face2node[iNode]];
			}
			face_pressure /= faceTopo->GetFaceNodeNum(iFace);
			force[0] += face_pressure * faceTopo->GetArea(iFace) * faceTopo->GetNormal(iFace)[0];
			force[1] += face_pressure * faceTopo->GetArea(iFace) * faceTopo->GetNormal(iFace)[1];
			force[2] += face_pressure * faceTopo->GetArea(iFace) * faceTopo->GetNormal(iFace)[2];
		}
		Log::info("Force: {}, {}, {}", force[0], force[1], force[2]);
	}



}