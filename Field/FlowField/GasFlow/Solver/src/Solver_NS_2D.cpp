#include"Solver_NS_2D.h"
#include"SpecialField.h"
namespace zaran
{
	void Solver_NS_2D::InitData()
	{
		// InitFieldShockReflection();
		// InitFieldIsentropicVortex();
		InitFieldNoFlow();
	}
	void Solver_NS_2D::InitFieldFarField()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		int init_type = para->GetInitFieldType();
		double prim_init[5];
		prim_init[0] = para->GetInflowDensity();
		prim_init[1] = para->GetInflowVelocityX();
		prim_init[2] = para->GetInflowVelocityY();
		prim_init[3] = para->GetInflowVelocityZ();
		prim_init[4] = para->GetInflowPressure();
		int total_node_num = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < total_node_num; ++iNode)
		{
			m_prim[0][iNode] = prim_init[0];
			m_prim[1][iNode] = prim_init[1];
			m_prim[2][iNode] = prim_init[2];
			m_prim[3][iNode] = prim_init[3];
			m_prim[4][iNode] = prim_init[4];
		}
		Prim2Cons();
	}
	void Solver_NS_2D::InitFieldNoFlow()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		int init_type = para->GetInitFieldType();
		double prim_init[5];
		prim_init[0] = para->GetInflowDensity();
		prim_init[1] = para->GetInflowVelocityX();
		prim_init[2] = para->GetInflowVelocityY();
		prim_init[3] = para->GetInflowVelocityZ();
		prim_init[4] = para->GetInflowPressure();
		int total_node_num = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < total_node_num; ++iNode)
		{
			m_prim[0][iNode] = 1 + 0.2 * coord[iNode].y();
			m_prim[1][iNode] = 3.0;
			m_prim[2][iNode] = 0;
			m_prim[3][iNode] = 0;
			m_prim[4][iNode] = prim_init[4];
		}
		Prim2Cons();
	}
	void Solver_NS_2D::InitFieldShockReflection()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		int init_type = para->GetInitFieldType();
		double prim_init[5];
		prim_init[0] = 6.4;
		prim_init[1] = 3.125;
		prim_init[2] = 0;
		prim_init[3] = 0;
		prim_init[4] = 18.5;
		int total_node_num = grid->GetTotalNodeNum();
		double x, y;
		for (int iNode = 0; iNode < total_node_num; ++iNode)
		{
			x = coord[iNode].x();
			y = coord[iNode].y();
			if (x <= 0.1)
			{
				m_prim[0][iNode] = prim_init[0];
				m_prim[1][iNode] = prim_init[1];
				m_prim[2][iNode] = prim_init[2];
				m_prim[3][iNode] = prim_init[3];
				m_prim[4][iNode] = prim_init[4];
			}
			else
			{
				m_prim[0][iNode] = 1.4;
				m_prim[1][iNode] = 0;
				m_prim[2][iNode] = 0;
				m_prim[3][iNode] = 0;
				m_prim[4][iNode] = 1.0;
			}
		}
		Prim2Cons();
	}
	void Solver_NS_2D::InitFieldIsentropicVortex()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		int total_node_num = grid->GetTotalNodeNum();
		double x, y;
		DVector prim;
		double beta = 5;
		for (int iNode = 0; iNode < total_node_num; ++iNode)
		{
			x = coord[iNode].x();
			y = coord[iNode].y();
			CalcIsentropicVortex(x, y, 5, prim);
			m_prim[0][iNode] = prim[0];
			m_prim[1][iNode] = prim[1];
			m_prim[2][iNode] = prim[2];
			m_prim[3][iNode] = prim[3];
			m_prim[4][iNode] = prim[4];
		}
		Prim2Cons();
	}
	void Solver_NS_2D::CalcMetric()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		auto& node_type = node_topo->GetType();
		auto& tempI = node_topo->GetTemplateI();
		auto& tempJ = node_topo->GetTemplateJ();
		int inner_node_num = grid->GetInnerNodeNum();
		CoordTrans coord_trans_coef;
		DVector3D x_right, x_left, y_right, y_left;
		for (size_t iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (node_type[iNode] != NodeType::inner && node_type[iNode] != NodeType::hole)
				continue;
			x_left = coord[tempI[iNode][0]];
			x_right = coord[tempI[iNode][2]];
			y_left = coord[tempJ[iNode][0]];
			y_right = coord[tempJ[iNode][2]];
			coord_trans_coef.CalcCoordTrans(int(grid->GetDimension()), x_right, x_left, y_right, y_left);
			if (coord_trans_coef.J() < 0)
			{
				tempJ[iNode] = IArray{ tempJ[iNode][2], tempJ[iNode][1], tempJ[iNode][0] };
				coord_trans_coef.CalcCoordTrans(int(grid->GetDimension()), x_right, x_left, y_left, y_right);
			}
			if (abs(coord_trans_coef.J()) < SMALL_NUMBER || isnan(abs(coord_trans_coef.J())) || isinf((coord_trans_coef.J())))
			{
				IArray currentTempI = { tempI[iNode][0], tempI[iNode][1], tempJ[iNode][0] };
				IArray currentTempJ = { tempI[iNode][2], tempJ[iNode][1], tempJ[iNode][0] };
				tempI[iNode] = currentTempI;
				tempJ[iNode] = currentTempJ;
				coord_trans_coef.CalcCoordTrans(int(grid->GetDimension()), x_left, y_left, x_right, y_right);
				if (coord_trans_coef.J() < 0)
				{
					tempJ[iNode] = IArray{ tempJ[iNode][2], tempJ[iNode][1], tempJ[iNode][0] };
					coord_trans_coef.CalcCoordTrans(int(grid->GetDimension()), x_left, y_left, y_right, x_right);
				}
			}

			// check coordinate
			if (coord_trans_coef.J() > 1e15 || coord_trans_coef.J() < 0)
			{
				Log::warn("Node {}: {},{},{}", iNode, coord[iNode].x(), coord[iNode].y(), coord[iNode].z());
				Log::info("x_left index={}: {},{},{}", tempI[iNode][0], x_left.x(), x_left.y(), x_left.z());
				Log::info("x_right index={}: {},{},{}", tempI[iNode][2], x_right.x(), x_right.y(), x_right.z());
				Log::info("y_left index={}: {},{},{}", tempJ[iNode][0], y_left.x(), y_left.y(), y_left.z());
				Log::info("y_right index={}: {},{},{}", tempJ[iNode][2], y_right.x(), y_right.y(), y_right.z());
			}
			m_metric[0][iNode] = coord_trans_coef.GetX()[0];
			m_metric[1][iNode] = coord_trans_coef.GetX()[1];
			m_metric[2][iNode] = coord_trans_coef.GetX()[2];
			m_metric[3][iNode] = coord_trans_coef.GetX()[3];
			m_metric[4][iNode] = coord_trans_coef.GetY()[0];
			m_metric[5][iNode] = coord_trans_coef.GetY()[1];
			m_metric[6][iNode] = coord_trans_coef.GetY()[2];
			m_metric[7][iNode] = coord_trans_coef.GetY()[3];
			m_metric[8][iNode] = coord_trans_coef.GetZ()[0];
			m_metric[9][iNode] = coord_trans_coef.GetZ()[1];
			m_metric[10][iNode] = coord_trans_coef.GetZ()[2];
			m_metric[11][iNode] = coord_trans_coef.GetZ()[3];
			m_metric[12][iNode] = coord_trans_coef.GetT()[0];
			m_metric[13][iNode] = coord_trans_coef.GetT()[1];
			m_metric[14][iNode] = coord_trans_coef.GetT()[2];
			m_metric[15][iNode] = coord_trans_coef.GetT()[3];
			m_metric[16][iNode] = coord_trans_coef.GetXi()[0];
			m_metric[17][iNode] = coord_trans_coef.GetXi()[1];
			m_metric[18][iNode] = coord_trans_coef.GetXi()[2];
			m_metric[19][iNode] = coord_trans_coef.GetXi()[3];
			m_metric[20][iNode] = coord_trans_coef.GetEta()[0];
			m_metric[21][iNode] = coord_trans_coef.GetEta()[1];
			m_metric[22][iNode] = coord_trans_coef.GetEta()[2];
			m_metric[23][iNode] = coord_trans_coef.GetEta()[3];
			m_metric[24][iNode] = coord_trans_coef.GetZeta()[0];
			m_metric[25][iNode] = coord_trans_coef.GetZeta()[1];
			m_metric[26][iNode] = coord_trans_coef.GetZeta()[2];
			m_metric[27][iNode] = coord_trans_coef.GetZeta()[3];
			m_metric[28][iNode] = coord_trans_coef.GetTau()[0];
			m_metric[29][iNode] = coord_trans_coef.GetTau()[1];
			m_metric[30][iNode] = coord_trans_coef.GetTau()[2];
			m_metric[31][iNode] = coord_trans_coef.GetTau()[3];
			m_metric[32][iNode] = coord_trans_coef.J();
		}
	}
	void Solver_NS_2D::CalcGradWLS()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		auto& coord = node_topo->GetCoordinate();
		auto& neighbor = node_topo->GetNeighborCloud();
		auto& prim = m_prim;
		auto& limiter_coef = m_limiter;
		auto& prim_grad_x = m_PrimGradX;
		auto& prim_grad_y = m_PrimGradY;
		auto& prim_grad_z = m_PrimGradZ;
		int inner_node_num = grid->GetInnerNodeNum();
		int bound_node_num = grid->GetBoundNodeNum();
		Matrix2d A;
		DVector2D b, grad;
		double omega = 0;
		double delta_value;
		double delta_x, delta_y;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (node_type[iNode] == NodeType::undefined)
				continue;
			auto& current_coord = coord[iNode];
			auto& current_neighbor = neighbor[iNode];
			for (size_t iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				A.setZero();
				b.setZero();
				for (size_t iNeib = 0; iNeib < current_neighbor.size(); ++iNeib)
				{
					omega = DistanceOfTwoPoints(coord[current_neighbor[iNeib]].data(), coord[iNode].data());
					if (abs(omega) < SMALL_NUMBER)
						continue;
					omega = 1.0 / omega;
					delta_value = m_prim[iVal][current_neighbor[iNeib]] - m_prim[iVal][iNode];
					delta_x = coord[current_neighbor[iNeib]].x() - coord[iNode].x();
					delta_y = coord[current_neighbor[iNeib]].y() - coord[iNode].y();
					A(0, 0) += omega * delta_x * delta_x;
					A(0, 1) += omega * delta_x * delta_y;
					A(1, 0) += omega * delta_y * delta_x;
					A(1, 1) += omega * delta_y * delta_y;
					b(0) += omega * delta_value * delta_x;
					b(1) += omega * delta_value * delta_y;
				}
				grad = A.inverse() * b;
				m_PrimGradX[iVal][iNode] = grad.x();
				m_PrimGradY[iVal][iNode] = grad.y();
			}
		}
	}
	void Solver_NS_2D::CalcTimeStepLocal()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		FlowSolverPara* para = GetPara();
		double cfl = para->GetCflNumber();
		int inner_node_num = grid->GetInnerNodeNum();
		double min_dt = LARGE_NUMBER;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{

			double gamma = 1.4;
			if (node_type[iNode] != NodeType::inner && node_type[iNode] != NodeType::hole)
				continue;
			double c = sqrt(gamma * m_prim[4][iNode] / m_prim[0][iNode]);
			double normXi = sqrt(m_metric[16][iNode] * m_metric[16][iNode] + m_metric[17][iNode] * m_metric[17][iNode] + m_metric[18][iNode] * m_metric[18][iNode]);
			double normEta = sqrt(m_metric[20][iNode] * m_metric[20][iNode] + m_metric[21][iNode] * m_metric[21][iNode] + m_metric[22][iNode] * m_metric[22][iNode]);
			double uXi = m_prim[1][iNode] * m_metric[16][iNode] + m_prim[2][iNode] * m_metric[17][iNode] + m_prim[3][iNode] * m_metric[18][iNode];
			double uEta = m_prim[1][iNode] * m_metric[20][iNode] + m_prim[2][iNode] * m_metric[21][iNode] + m_prim[3][iNode] * m_metric[22][iNode];
			double lamda = abs(uXi) + abs(uEta) + c * (normXi + normEta);
			m_dt[iNode] = cfl / lamda;
			min_dt = Min(min_dt, m_dt[iNode]);
		}
		GlobalData::Update("dt", min_dt);
	}

	void Solver_NS_2D::InviscidFlux()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		auto& node_coord = node_topo->GetCoordinate();
		auto& template_i = node_topo->GetTemplateI();
		auto& template_j = node_topo->GetTemplateJ();
		auto& res = m_residual;
		RiemannSolverPara riemann_para;
		riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
		for (size_t iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (node_type[iNode] != NodeType::inner)
				continue;
			auto& jacobi = m_metric[32][iNode];
			// i direction
			riemann_para.norm(0) = m_metric[16][iNode];
			riemann_para.norm(1) = m_metric[17][iNode];
			riemann_para.norm(2) = m_metric[18][iNode];
			riemann_para.nt = m_metric[19][iNode];
			MidPointReconstruct(template_i[iNode][1], template_i[iNode][2], &riemann_para.prim_left(0), &riemann_para.prim_right(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				m_residual[iVal][iNode] += riemann_para.flux[iVal] / jacobi;
			MidPointReconstruct(template_i[iNode][0], template_i[iNode][1], &riemann_para.prim_left(0), &riemann_para.prim_right(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				m_residual[iVal][iNode] -= riemann_para.flux[iVal] / jacobi;
			// j direction
			riemann_para.norm(0) = m_metric[20][iNode];
			riemann_para.norm(1) = m_metric[21][iNode];
			riemann_para.norm(2) = m_metric[22][iNode];
			riemann_para.nt = m_metric[23][iNode];
			MidPointReconstruct(template_j[iNode][1], template_j[iNode][2], &riemann_para.prim_left(0), &riemann_para.prim_right(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				m_residual[iVal][iNode] += riemann_para.flux[iVal] / jacobi;
			MidPointReconstruct(template_j[iNode][0], template_j[iNode][1], &riemann_para.prim_left(0), &riemann_para.prim_right(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				m_residual[iVal][iNode] -= riemann_para.flux[iVal] / jacobi;
		}
	}

	void Solver_NS_2D::SolveHoleNode()
	{
		string hole_solver_type = GlobalData::GetString("HoleSolverType");
		if (hole_solver_type == "FNFDM")
			SolveHoleNodeFNFDM();
		else if (hole_solver_type == "IDW")
			SolveHoleNodeIDW();
		else
			Log::error("HoleSolverType: {} is not defined!", hole_solver_type);
	}
	void Solver_NS_2D::SolveHoleNodeFNFDM()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		auto& template_i = node_topo->GetTemplateI();
		auto& template_j = node_topo->GetTemplateJ();
		auto& neibor = node_topo->GetNeighborCloud();
		auto& coord = node_topo->GetCoordinate();
		int inner_node_num = grid->GetInnerNodeNum();
		int bound_node_num = grid->GetBoundNodeNum();
		DVector2D r;
		RiemannSolverPara riemann_para;
		riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (node_type[iNode] != NodeType::hole)
				continue;
			auto& jacobi = m_metric[32][iNode];
			// i direction
			riemann_para.norm(0) = m_metric[16][iNode];
			riemann_para.norm(1) = m_metric[17][iNode];
			riemann_para.norm(2) = 0;
			riemann_para.nt = m_metric[19][iNode];
			IArray tempI = template_i[iNode];
			r[0] = coord[tempI[2]][0] - coord[iNode][0];
			r[1] = coord[tempI[2]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.prim_left(iVal) = m_prim[iVal][iNode] + 0.5 * m_limiter[iVal][iNode] * (r[0] * m_PrimGradX[iVal][iNode] + r[1] * m_PrimGradY[iVal][iNode]);
				riemann_para.prim_right(iVal) = m_prim[iVal][tempI[2]] - 0.5 * m_limiter[iVal][tempI[2]] * (r[0] * m_PrimGradX[iVal][tempI[2]] + r[1] * m_PrimGradY[iVal][tempI[2]]);
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				m_residual[iVal][iNode] += riemann_para.flux[iVal] / jacobi;
			r[0] = coord[tempI[0]][0] - coord[iNode][0];
			r[1] = coord[tempI[0]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.prim_left(iVal) = m_prim[iVal][tempI[0]] - 0.5 * m_limiter[iVal][tempI[0]] * (r[0] * m_PrimGradX[iVal][tempI[0]] + r[1] * m_PrimGradY[iVal][tempI[0]]);
				riemann_para.prim_right(iVal) = m_prim[iVal][iNode] + 0.5 * m_limiter[iVal][iNode] * (r[0] * m_PrimGradX[iVal][iNode] + r[1] * m_PrimGradY[iVal][iNode]);
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				m_residual[iVal][iNode] -= riemann_para.flux[iVal] / jacobi;

			// j direction
			riemann_para.norm(0) = m_metric[20][iNode];
			riemann_para.norm(1) = m_metric[21][iNode];
			riemann_para.norm(2) = 0;
			riemann_para.nt = m_metric[23][iNode];
			IArray tempJ = template_j[iNode];
			r[0] = coord[tempJ[2]][0] - coord[iNode][0];
			r[1] = coord[tempJ[2]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.prim_left(iVal) = m_prim[iVal][iNode] + 0.5 * m_limiter[iVal][iNode] * (r[0] * m_PrimGradX[iVal][iNode] + r[1] * m_PrimGradY[iVal][iNode]);
				riemann_para.prim_right(iVal) = m_prim[iVal][tempJ[2]] - 0.5 * m_limiter[iVal][tempJ[2]] * (r[0] * m_PrimGradX[iVal][tempJ[2]] + r[1] * m_PrimGradY[iVal][tempJ[2]]);
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				m_residual[iVal][iNode] += riemann_para.flux[iVal] / jacobi;
			r[0] = coord[tempJ[0]][0] - coord[iNode][0];
			r[1] = coord[tempJ[0]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.prim_left(iVal) = m_prim[iVal][tempJ[0]] - 0.5 * m_limiter[iVal][tempJ[0]] * (r[0] * m_PrimGradX[iVal][tempJ[0]] + r[1] * m_PrimGradY[iVal][tempJ[0]]);
				riemann_para.prim_right(iVal) = m_prim[iVal][iNode] + 0.5 * m_limiter[iVal][iNode] * (r[0] * m_PrimGradX[iVal][iNode] + r[1] * m_PrimGradY[iVal][iNode]);
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				m_residual[iVal][iNode] -= riemann_para.flux[iVal] / jacobi;
				if (isnan(m_residual[iVal][iNode]) || isinf(m_residual[iVal][iNode]))
					Log::error("inode={},NAN in Residual!", iNode);
				m_cons[iVal][iNode] -= m_residual[iVal][iNode] * m_dt[iNode] * m_metric[32][iNode];
				m_residual[iVal][iNode] = 0;
			}
		}
	}
	void Solver_NS_2D::SolveHoleNodeIDW()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		auto& neibor = node_topo->GetNeighborCloud();
		auto& coord = node_topo->GetCoordinate();
		auto& template_i = node_topo->GetTemplateI();
		auto& template_j = node_topo->GetTemplateJ();
		Array<int> neiborNode;
		DArray weight, distance;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (node_type[iNode] != NodeType::hole)
				continue;
			auto& tempI = template_i[iNode];
			auto& tempJ = template_j[iNode];
			neiborNode = { tempI[0],tempI[2],tempJ[0],tempJ[2] };
			weight.resize(neiborNode.size());
			distance.resize(neiborNode.size());
			double sum = 0;
			for (int iNeib = 0; iNeib < neiborNode.size(); ++iNeib)
			{
				distance[iNeib] = DistanceOfTwoPoints(coord[iNode].data(), coord[neiborNode[iNeib]].data());
				sum += 1.0 / distance[iNeib];
			}
			for (int iNeib = 0; iNeib < neiborNode.size(); ++iNeib)
				weight[iNeib] = 1.0 / distance[iNeib] / sum;
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				m_prim[iVal][iNode] = 0;
				m_residual[iVal][iNode] = 0;
				for (int iNeib = 0; iNeib < neiborNode.size(); ++iNeib)
				{
					m_prim[iVal][iNode] += m_prim[iVal][neiborNode[iNeib]] * weight[iNeib];
				}
			}
			Prim2Cons(m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode], m_cons[0][iNode], m_cons[1][iNode], m_cons[2][iNode], m_cons[3][iNode], m_cons[4][iNode]);
		}
	}
	void Solver_NS_2D::BoundaryCondition()
	{
		NSSolver::BoundaryCondition();
		SolveHoleNode();
		SolveUserDefinedBoundary();
	}
	void Solver_NS_2D::SolveUserDefinedBoundary()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		auto& neibor = node_topo->GetNeighborCloud();
		auto& coord = node_topo->GetCoordinate();
		auto& template_i = node_topo->GetTemplateI();
		auto& template_j = node_topo->GetTemplateJ();
		double x, y;
		DVector prim_ideal(5);
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (node_type[iNode] != NodeType::userDefined)
				continue;
			x = coord[iNode].x();
			y = coord[iNode].y();
			CalcIsentropicVortex(x, y, 5.0, prim_ideal);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				m_prim[iVal][iNode] = prim_ideal(iVal);
				m_residual[iVal][iNode] = 0;
			}
			Prim2Cons(m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode], m_cons[0][iNode], m_cons[1][iNode], m_cons[2][iNode], m_cons[3][iNode], m_cons[4][iNode]);
		}
	}
	void Solver_NS_2D::ViscousFlux()
	{

	}

	void Solver_NS_2D::SourceFlux()
	{
	}

}