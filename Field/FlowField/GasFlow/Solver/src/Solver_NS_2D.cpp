#include"Solver_NS_2D.h"
#include"SpecialField.h"
namespace zaran
{
	void Solver_NS_2D::InitField()
	{
		// InitFieldShockReflection();
		// InitFieldIsentropicVortex();
		InitFieldNoFlow();
	}
	void Solver_NS_2D::InitFieldFarFlow()
	{
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		const InitFieldType& init_type = para->GetInitFieldType();
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
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		const InitFieldType& init_type = para->GetInitFieldType();
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
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		FlowSolverPara* para = GetPara();
		const InitFieldType& init_type = para->GetInitFieldType();
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
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
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
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
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
			if (coord_trans_coef.Jacobian() < 0)
			{
				tempJ[iNode] = IArray{ tempJ[iNode][2], tempJ[iNode][1], tempJ[iNode][0] };
				coord_trans_coef.CalcCoordTrans(int(grid->GetDimension()), x_right, x_left, y_left, y_right);
			}
			if (abs(coord_trans_coef.Jacobian()) < SMALL_NUMBER || isnan(abs(coord_trans_coef.Jacobian())) || isinf((coord_trans_coef.Jacobian())))
			{
				IArray currentTempI = { tempI[iNode][0], tempI[iNode][1], tempJ[iNode][0] };
				IArray currentTempJ = { tempI[iNode][2], tempJ[iNode][1], tempJ[iNode][0] };
				tempI[iNode] = currentTempI;
				tempJ[iNode] = currentTempJ;
				coord_trans_coef.CalcCoordTrans(int(grid->GetDimension()), x_left, y_left, x_right, y_right);
				if (coord_trans_coef.Jacobian() < 0)
				{
					tempJ[iNode] = IArray{ tempJ[iNode][2], tempJ[iNode][1], tempJ[iNode][0] };
					coord_trans_coef.CalcCoordTrans(int(grid->GetDimension()), x_left, y_left, y_right, x_right);
				}
			}

			// check coordinate
			if (coord_trans_coef.Jacobian() > 1e15 || coord_trans_coef.Jacobian() < 0)
			{
				Log::warn("Node {}: {},{},{}", iNode, coord[iNode].x(), coord[iNode].y(), coord[iNode].z());
				Log::info("x_left index={}: {},{},{}", tempI[iNode][0], x_left.x(), x_left.y(), x_left.z());
				Log::info("x_right index={}: {},{},{}", tempI[iNode][2], x_right.x(), x_right.y(), x_right.z());
				Log::info("y_left index={}: {},{},{}", tempJ[iNode][0], y_left.x(), y_left.y(), y_left.z());
				Log::info("y_right index={}: {},{},{}", tempJ[iNode][2], y_right.x(), y_right.y(), y_right.z());
			}
			GetMetricXi(iNode)[0] = coord_trans_coef.GetXi()[0];
			GetMetricXi(iNode)[1] = coord_trans_coef.GetXi()[1];
			GetMetricXi(iNode)[2] = coord_trans_coef.GetXi()[2];
			GetMetricXi(iNode)[3] = coord_trans_coef.GetXi()[3];
			GetMetricEta(iNode)[0] = coord_trans_coef.GetEta()[0];
			GetMetricEta(iNode)[1] = coord_trans_coef.GetEta()[1];
			GetMetricEta(iNode)[2] = coord_trans_coef.GetEta()[2];
			GetMetricEta(iNode)[3] = coord_trans_coef.GetEta()[3];
			GetMetricZeta(iNode)[0] = coord_trans_coef.GetZeta()[0];
			GetMetricZeta(iNode)[1] = coord_trans_coef.GetZeta()[1];
			GetMetricZeta(iNode)[2] = coord_trans_coef.GetZeta()[2];
			GetMetricZeta(iNode)[3] = coord_trans_coef.GetZeta()[3];
			GetMetricTau(iNode)[0] = coord_trans_coef.GetTau()[0];
			GetMetricTau(iNode)[1] = coord_trans_coef.GetTau()[1];
			GetMetricTau(iNode)[2] = coord_trans_coef.GetTau()[2];
			GetMetricTau(iNode)[3] = coord_trans_coef.GetTau()[3];
			GetMetricJacob(iNode) = coord_trans_coef.Jacobian();
		}
	}
	void Solver_NS_2D::CalcGradWLS()
	{
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		auto& coord = node_topo->GetCoordinate();
		auto& neighbor = node_topo->GetNeighborCloud();
		auto& prim = m_prim;
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
				GetPrimGrad(iNode, iVal, 0) = grad.x();
				GetPrimGrad(iNode, iVal, 1) = grad.y();
			}
		}
	}
	void Solver_NS_2D::CalcTimeStepLocal()
	{
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
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
			double normXi = sqrt(GetMetricXi(iNode)[0] * GetMetricXi(iNode)[0] + GetMetricXi(iNode)[1] * GetMetricXi(iNode)[1] + GetMetricXi(iNode)[2] * GetMetricXi(iNode)[2]);
			double normEta = sqrt(GetMetricEta(iNode)[0] * GetMetricEta(iNode)[0] + GetMetricEta(iNode)[1] * GetMetricEta(iNode)[1] + GetMetricEta(iNode)[2] * GetMetricEta(iNode)[2]);
			double uXi = m_prim[1][iNode] * GetMetricXi(iNode)[0] + m_prim[2][iNode] * GetMetricXi(iNode)[1] + m_prim[3][iNode] * GetMetricXi(iNode)[2];
			double uEta = m_prim[1][iNode] * GetMetricEta(iNode)[0] + m_prim[2][iNode] * GetMetricEta(iNode)[1] + m_prim[3][iNode] * GetMetricEta(iNode)[2];
			double lamda = abs(uXi) + abs(uEta) + c * (normXi + normEta);
			m_dt[iNode] = cfl / lamda;
			min_dt = Min(min_dt, m_dt[iNode]);
		}
		GlobalData::Update("dt", min_dt);
	}

	void Solver_NS_2D::ConvectiveResidual()
	{
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		auto& node_coord = node_topo->GetCoordinate();
		auto& template_i = node_topo->GetTemplateI();
		auto& template_j = node_topo->GetTemplateJ();
		RiemannSolverPara riemann_para;
		riemann_para.gamma_left = riemann_para.gamma_right = 1.4;
		for (size_t iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (node_type[iNode] != NodeType::inner)
				continue;
			auto& jacobi = GetMetricJacob(iNode);
			// i direction
			riemann_para.norm(0) = GetMetricXi(iNode)[0];
			riemann_para.norm(1) = GetMetricXi(iNode)[1];
			riemann_para.norm(2) = GetMetricXi(iNode)[2];
			riemann_para.nt = GetMetricXi(iNode)[3];
			MidPointReconstruct(template_i[iNode][1], template_i[iNode][2], &riemann_para.prim_left(0), &riemann_para.prim_right(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
			MidPointReconstruct(template_i[iNode][0], template_i[iNode][1], &riemann_para.prim_left(0), &riemann_para.prim_right(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;
			// j direction
			riemann_para.norm(0) = GetMetricEta(iNode)[0];
			riemann_para.norm(1) = GetMetricEta(iNode)[1];
			riemann_para.norm(2) = GetMetricEta(iNode)[2];
			riemann_para.nt = GetMetricEta(iNode)[3];
			MidPointReconstruct(template_j[iNode][1], template_j[iNode][2], &riemann_para.prim_left(0), &riemann_para.prim_right(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
			MidPointReconstruct(template_j[iNode][0], template_j[iNode][1], &riemann_para.prim_left(0), &riemann_para.prim_right(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;
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
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
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
			auto& jacobi = GetMetricJacob(iNode);
			// i direction
			riemann_para.norm(0) = GetMetricXi(iNode)[0];
			riemann_para.norm(1) = GetMetricXi(iNode)[1];
			riemann_para.norm(2) = 0;
			riemann_para.nt = GetMetricXi(iNode)[3];
			IArray tempI = template_i[iNode];
			r[0] = coord[tempI[2]][0] - coord[iNode][0];
			r[1] = coord[tempI[2]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.prim_left(iVal) = m_prim[iVal][iNode] + 0.5 * GetLimiter(iNode, iVal) * (r[0] * GetPrimGrad(iNode, iVal, 0) + r[1] * GetPrimGrad(iNode, iVal, 1));
				riemann_para.prim_right(iVal) = m_prim[iVal][tempI[2]] - 0.5 * GetLimiter(tempI[2], iVal) * (r[0] * GetPrimGrad(tempI[2], iVal, 0) + r[1] * GetPrimGrad(tempI[2], iVal, 1));
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
			r[0] = coord[tempI[0]][0] - coord[iNode][0];
			r[1] = coord[tempI[0]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.prim_left(iVal) = m_prim[iVal][tempI[0]] - 0.5 * GetLimiter(tempI[0], iVal) * (r[0] * GetPrimGrad(tempI[0], iVal, 0) + r[1] * GetPrimGrad(tempI[0], iVal, 1));
				riemann_para.prim_right(iVal) = m_prim[iVal][iNode] + 0.5 * GetLimiter(iNode, iVal) * (r[0] * GetPrimGrad(iNode, iVal, 0) + r[1] * GetPrimGrad(iNode, iVal, 1));
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;

			// j direction
			riemann_para.norm(0) = GetMetricEta(iNode)[0];
			riemann_para.norm(1) = GetMetricEta(iNode)[1];
			riemann_para.norm(2) = 0;
			riemann_para.nt = GetMetricEta(iNode)[3];
			IArray tempJ = template_j[iNode];
			r[0] = coord[tempJ[2]][0] - coord[iNode][0];
			r[1] = coord[tempJ[2]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.prim_left(iVal) = m_prim[iVal][iNode] + 0.5 * GetLimiter(iNode, iVal) * (r[0] * GetPrimGrad(iNode, iVal, 0) + r[1] * GetPrimGrad(iNode, iVal, 1));
				riemann_para.prim_right(iVal) = m_prim[iVal][tempJ[2]] - 0.5 * GetLimiter(tempJ[2], iVal) * (r[0] * GetPrimGrad(tempJ[2], iVal, 0) + r[1] * GetPrimGrad(tempJ[2], iVal, 1));
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				GetResidual(iNode, iVal) += riemann_para.flux[iVal] / jacobi;
			r[0] = coord[tempJ[0]][0] - coord[iNode][0];
			r[1] = coord[tempJ[0]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.prim_left(iVal) = m_prim[iVal][tempJ[0]] - 0.5 * GetLimiter(tempJ[0], iVal) * (r[0] * GetPrimGrad(tempJ[0], iVal, 0) + r[1] * GetPrimGrad(tempJ[0], iVal, 1));
				riemann_para.prim_right(iVal) = m_prim[iVal][iNode] + 0.5 * GetLimiter(iNode, iVal) * (r[0] * GetPrimGrad(iNode, iVal, 0) + r[1] * GetPrimGrad(iNode, iVal, 1));
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				GetResidual(iNode, iVal) -= riemann_para.flux[iVal] / jacobi;
				if (isnan(GetResidual(iNode, iVal)) || isinf(GetResidual(iNode, iVal)))
					Log::error("inode={},NAN in Residual!", iNode);
				GetCons(iNode, iVal) -= GetResidual(iNode, iVal) * m_dt[iNode] * GetMetricJacob(iNode);
				GetResidual(iNode, iVal) = 0;
			}
		}
	}
	void Solver_NS_2D::SolveHoleNodeIDW()
	{
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
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
				GetResidual(iNode, iVal) = 0;
				for (int iNeib = 0; iNeib < neiborNode.size(); ++iNeib)
				{
					m_prim[iVal][iNode] += m_prim[iVal][neiborNode[iNeib]] * weight[iNeib];
				}
			}
			Prim2Cons(m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode], GetCons(iNode, 0), GetCons(iNode, 1), GetCons(iNode, 2), GetCons(iNode, 3), GetCons(iNode, 4));
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
		Grid* grid = GetGrid();
		NodeTopo* node_topo = grid->GetNodeTopo();
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
				GetResidual(iNode, iVal) = 0;
			}
			Prim2Cons(m_prim[0][iNode], m_prim[1][iNode], m_prim[2][iNode], m_prim[3][iNode], m_prim[4][iNode], GetCons(iNode, 0), GetCons(iNode, 1), GetCons(iNode, 2), GetCons(iNode, 3), GetCons(iNode, 4));
		}
	}
	void Solver_NS_2D::ViscousResidual()
	{

	}

	void Solver_NS_2D::SourceTermResidual()
	{
	}

}