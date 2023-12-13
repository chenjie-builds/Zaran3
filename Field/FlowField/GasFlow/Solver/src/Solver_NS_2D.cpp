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
	void Solver_NS_2D::InitFieldFarField()
	{
		GridPtr grid = GetGrid();
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
		auto& node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		FlowSolverParaPtr para = GetPara();
		int init_type = para->GetInitFieldType();
		DVector prim_init = para->GetPrimitiveInflow();
		int total_node_num = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < total_node_num; ++iNode)
		{
			rho[iNode] = prim_init[0];
			u[iNode] = prim_init[1];
			v[iNode] = prim_init[2];
			w[iNode] = prim_init[3];
			p[iNode] = prim_init[4];
		}
		Primitive2Conservative();
	}
	void Solver_NS_2D::InitFieldNoFlow()
	{
		GridPtr grid = GetGrid();
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
		auto& node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		FlowSolverParaPtr para = GetPara();
		int init_type = para->GetInitFieldType();
		DVector prim_init = para->GetPrimitiveInflow();
		int total_node_num = grid->GetTotalNodeNum();
		for (int iNode = 0; iNode < total_node_num; ++iNode)
		{
			rho[iNode] = 1 + 0.2 * coord[iNode].y();
			u[iNode] = 3.0;
			v[iNode] = 0;
			w[iNode] = 0;
			p[iNode] = prim_init[4];
		}
		Primitive2Conservative();
	}
	void Solver_NS_2D::InitFieldShockReflection()
	{
		GridPtr grid = GetGrid();
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
		auto& node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		FlowSolverParaPtr para = GetPara();
		int init_type = para->GetInitFieldType();
		DVector prim_init = para->GetPrimitiveInflow();
		prim_init[0] = 6.4;
		prim_init[1] = 3.125;
		prim_init[2] = 0;
		prim_init[3] = 0;
		prim_init[4] = 18.5;
		int total_node_num = grid->GetTotalNodeNum();
		double x, y;
		para->SetPrimitiveInflow(prim_init);
		for (int iNode = 0; iNode < total_node_num; ++iNode)
		{
			x = coord[iNode].x();
			y = coord[iNode].y();
			if (x <= 0.1)
			{
				rho[iNode] = prim_init[0];
				u[iNode] = prim_init[1];
				v[iNode] = prim_init[2];
				w[iNode] = prim_init[3];
				p[iNode] = prim_init[4];
			}
			else
			{
				rho[iNode] = 1.4;
				u[iNode] = 0;
				v[iNode] = 0;
				w[iNode] = 0;
				p[iNode] = 1.0;
			}
		}
		Primitive2Conservative();
	}
	void Solver_NS_2D::InitFieldIsentropicVortex()
	{
		GridPtr grid = GetGrid();
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
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
			rho[iNode] = prim[0];
			u[iNode] = prim[1];
			v[iNode] = prim[2];
			w[iNode] = prim[3];
			p[iNode] = prim[4];
		}
		Primitive2Conservative();
	}
	void Solver_NS_2D::ComputeCoordTrans()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& coord = node_topo->GetCoordinate();
		auto& node_type = node_topo->GetType();
		auto& tempI = node_topo->GetTemplateI();
		auto& tempJ = node_topo->GetTemplateJ();
		auto& coordTransCoef = m_CoordTrans;
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
				ZaranLog::warn("Node {}: {},{},{}", iNode, coord[iNode].x(), coord[iNode].y(), coord[iNode].z());
				ZaranLog::info("x_left index={}: {},{},{}", tempI[iNode][0], x_left.x(), x_left.y(), x_left.z());
				ZaranLog::info("x_right index={}: {},{},{}", tempI[iNode][2], x_right.x(), x_right.y(), x_right.z());
				ZaranLog::info("y_left index={}: {},{},{}", tempJ[iNode][0], y_left.x(), y_left.y(), y_left.z());
				ZaranLog::info("y_right index={}: {},{},{}", tempJ[iNode][2], y_right.x(), y_right.y(), y_right.z());
			}

			(*coordTransCoef[0])[iNode] = coord_trans_coef.GetX()[0];
			(*coordTransCoef[1])[iNode] = coord_trans_coef.GetX()[1];
			(*coordTransCoef[2])[iNode] = coord_trans_coef.GetX()[2];
			(*coordTransCoef[3])[iNode] = coord_trans_coef.GetX()[3];
			(*coordTransCoef[4])[iNode] = coord_trans_coef.GetY()[0];
			(*coordTransCoef[5])[iNode] = coord_trans_coef.GetY()[1];
			(*coordTransCoef[6])[iNode] = coord_trans_coef.GetY()[2];
			(*coordTransCoef[7])[iNode] = coord_trans_coef.GetY()[3];
			(*coordTransCoef[8])[iNode] = coord_trans_coef.GetZ()[0];
			(*coordTransCoef[9])[iNode] = coord_trans_coef.GetZ()[1];
			(*coordTransCoef[10])[iNode] = coord_trans_coef.GetZ()[2];
			(*coordTransCoef[11])[iNode] = coord_trans_coef.GetZ()[3];
			(*coordTransCoef[12])[iNode] = coord_trans_coef.GetT()[0];
			(*coordTransCoef[13])[iNode] = coord_trans_coef.GetT()[1];
			(*coordTransCoef[14])[iNode] = coord_trans_coef.GetT()[2];
			(*coordTransCoef[15])[iNode] = coord_trans_coef.GetT()[3];
			(*coordTransCoef[16])[iNode] = coord_trans_coef.GetXi()[0];
			(*coordTransCoef[17])[iNode] = coord_trans_coef.GetXi()[1];
			(*coordTransCoef[18])[iNode] = coord_trans_coef.GetXi()[2];
			(*coordTransCoef[19])[iNode] = coord_trans_coef.GetXi()[3];
			(*coordTransCoef[20])[iNode] = coord_trans_coef.GetEta()[0];
			(*coordTransCoef[21])[iNode] = coord_trans_coef.GetEta()[1];
			(*coordTransCoef[22])[iNode] = coord_trans_coef.GetEta()[2];
			(*coordTransCoef[23])[iNode] = coord_trans_coef.GetEta()[3];
			(*coordTransCoef[24])[iNode] = coord_trans_coef.GetZeta()[0];
			(*coordTransCoef[25])[iNode] = coord_trans_coef.GetZeta()[1];
			(*coordTransCoef[26])[iNode] = coord_trans_coef.GetZeta()[2];
			(*coordTransCoef[27])[iNode] = coord_trans_coef.GetZeta()[3];
			(*coordTransCoef[28])[iNode] = coord_trans_coef.GetTau()[0];
			(*coordTransCoef[29])[iNode] = coord_trans_coef.GetTau()[1];
			(*coordTransCoef[30])[iNode] = coord_trans_coef.GetTau()[2];
			(*coordTransCoef[31])[iNode] = coord_trans_coef.GetTau()[3];
			(*coordTransCoef[32])[iNode] = coord_trans_coef.J();
		}
	}
	void Solver_NS_2D::ComputeGradientWLS()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		auto& coord = node_topo->GetCoordinate();
		auto& neighbor = node_topo->GetNeighborCloud();
		auto& prim = m_Primitive;
		auto& limiter_coef = m_LimiterCoef;
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
					delta_value = (*prim[iVal])[current_neighbor[iNeib]] - (*prim[iVal])[iNode];
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
				(*prim_grad_x[iVal])[iNode] = grad.x();
				(*prim_grad_y[iVal])[iNode] = grad.y();
			}
		}
	}
	void Solver_NS_2D::ComputeTimeStepLocal()
	{
		GridPtr grid = GetGrid();
		auto& node_topo = grid->GetNodeTopo();
		auto& node_type = node_topo->GetType();
		FlowSolverParaPtr para = GetPara();
		auto& rho = *m_Primitive[0];
		auto& u = *m_Primitive[1];
		auto& v = *m_Primitive[2];
		auto& w = *m_Primitive[3];
		auto& p = *m_Primitive[4];
		auto& dt = *m_TimeStep;
		auto& coord_trans_coef = m_CoordTrans;
		double cfl = para->GetCflNumber();
		int inner_node_num = grid->GetInnerNodeNum();
		double min_dt = LARGE_NUMBER;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{

			double gamma = 1.4;
			if (node_type[iNode] != NodeType::inner && node_type[iNode] != NodeType::hole)
				continue;
			double c = sqrt(gamma * p[iNode] / rho[iNode]);
			double normXi = sqrt((*coord_trans_coef[16])[iNode] * (*coord_trans_coef[16])[iNode] + (*coord_trans_coef[17])[iNode] * (*coord_trans_coef[17])[iNode] + (*coord_trans_coef[18])[iNode] * (*coord_trans_coef[19])[iNode]);
			double normEta = sqrt((*coord_trans_coef[20])[iNode] * (*coord_trans_coef[20])[iNode] + (*coord_trans_coef[21])[iNode] * (*coord_trans_coef[21])[iNode] + (*coord_trans_coef[22])[iNode] * (*coord_trans_coef[22])[iNode]);
			double uXi = u[iNode] * (*coord_trans_coef[16])[iNode] + v[iNode] * (*coord_trans_coef[17])[iNode] + w[iNode] * (*coord_trans_coef[18])[iNode] + (*coord_trans_coef[19])[iNode];
			double uEta = u[iNode] * (*coord_trans_coef[20])[iNode] + v[iNode] * (*coord_trans_coef[21])[iNode] + w[iNode] * (*coord_trans_coef[22])[iNode] + (*coord_trans_coef[23])[iNode];
			double lamda = abs(uXi) + abs(uEta) + c * (normXi + normEta);
			dt[iNode] = cfl / lamda;
			min_dt = Min(min_dt, dt[iNode]);
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
		auto& res = m_Residual;
		RiemannSolverPara riemann_para;
		riemann_para.gammaL = riemann_para.gammaR = 1.4;
		for (size_t iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (node_type[iNode] != NodeType::inner)
				continue;
			auto& jacobi = (*m_CoordTrans[32])[iNode];
			// i direction
			riemann_para.norm(0) = (*m_CoordTrans[16])[iNode];
			riemann_para.norm(1) = (*m_CoordTrans[17])[iNode];
			riemann_para.norm(2) = (*m_CoordTrans[18])[iNode];
			riemann_para.nt = (*m_CoordTrans[19])[iNode];
			MidPointReconstruct(template_i[iNode][1], template_i[iNode][2], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
			MidPointReconstruct(template_i[iNode][0], template_i[iNode][1], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
			// j direction
			riemann_para.norm(0) = (*m_CoordTrans[20])[iNode];
			riemann_para.norm(1) = (*m_CoordTrans[21])[iNode];
			riemann_para.norm(2) = (*m_CoordTrans[22])[iNode];
			riemann_para.nt = (*m_CoordTrans[23])[iNode];
			MidPointReconstruct(template_j[iNode][1], template_j[iNode][2], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
			MidPointReconstruct(template_j[iNode][0], template_j[iNode][1], &riemann_para.primL(0), &riemann_para.primR(0));
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
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
			ZaranLog::error("HoleSolverType: {} is not defined!", hole_solver_type);
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
		auto& prim = m_Primitive;
		auto& cons = m_Conservative;
		auto& prim_grad_x = m_PrimGradX;
		auto& prim_grad_y = m_PrimGradY;
		auto& limiter_coef = m_LimiterCoef;
		auto& res = m_Residual;
		auto& dt = *m_TimeStep;
		auto& coord_trans_coef = m_CoordTrans;
		int inner_node_num = grid->GetInnerNodeNum();
		int bound_node_num = grid->GetBoundNodeNum();
		DVector2D r;
		RiemannSolverPara riemann_para;
		riemann_para.gammaL = riemann_para.gammaR = 1.4;
		for (int iNode = 0; iNode < grid->GetTotalNodeNum(); ++iNode)
		{
			if (node_type[iNode] != NodeType::hole)
				continue;
			auto& jacobi = (*coord_trans_coef[32])[iNode];
			// i direction
			riemann_para.norm(0) = (*coord_trans_coef[16])[iNode];
			riemann_para.norm(1) = (*coord_trans_coef[17])[iNode];
			riemann_para.norm(2) = 0;
			riemann_para.nt = (*coord_trans_coef[19])[iNode];
			IArray tempI = template_i[iNode];
			r[0] = coord[tempI[2]][0] - coord[iNode][0];
			r[1] = coord[tempI[2]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiter_coef[iVal])[iNode] * (r[0] * (*prim_grad_x[iVal])[iNode] + r[1] * (*prim_grad_y[iVal])[iNode]);
				riemann_para.primR(iVal) = (*prim[iVal])[tempI[2]] - 0.5 * (*limiter_coef[iVal])[tempI[2]] * (r[0] * (*prim_grad_x[iVal])[tempI[2]] + r[1] * (*prim_grad_y[iVal])[tempI[2]]);
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
			r[0] = coord[tempI[0]][0] - coord[iNode][0];
			r[1] = coord[tempI[0]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.primL(iVal) = (*prim[iVal])[tempI[0]] - 0.5 * (*limiter_coef[iVal])[tempI[0]] * (r[0] * (*prim_grad_x[iVal])[tempI[0]] + r[1] * (*prim_grad_y[iVal])[tempI[0]]);
				riemann_para.primR(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiter_coef[iVal])[iNode] * (r[0] * (*prim_grad_x[iVal])[iNode] + r[1] * (*prim_grad_y[iVal])[iNode]);
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;

			// j direction
			riemann_para.norm(0) = (*coord_trans_coef[20])[iNode];
			riemann_para.norm(1) = (*coord_trans_coef[21])[iNode];
			riemann_para.norm(2) = 0;
			riemann_para.nt = (*coord_trans_coef[23])[iNode];
			IArray tempJ = template_j[iNode];
			r[0] = coord[tempJ[2]][0] - coord[iNode][0];
			r[1] = coord[tempJ[2]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.primL(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiter_coef[iVal])[iNode] * (r[0] * (*prim_grad_x[iVal])[iNode] + r[1] * (*prim_grad_y[iVal])[iNode]);
				riemann_para.primR(iVal) = (*prim[iVal])[tempJ[2]] - 0.5 * (*limiter_coef[iVal])[tempJ[2]] * (r[0] * (*prim_grad_x[iVal])[tempJ[2]] + r[1] * (*prim_grad_y[iVal])[tempJ[2]]);
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
				(*res[iVal])[iNode] += riemann_para.flux[iVal] / jacobi;
			r[0] = coord[tempJ[0]][0] - coord[iNode][0];
			r[1] = coord[tempJ[0]][1] - coord[iNode][1];
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				riemann_para.primL(iVal) = (*prim[iVal])[tempJ[0]] - 0.5 * (*limiter_coef[iVal])[tempJ[0]] * (r[0] * (*prim_grad_x[iVal])[tempJ[0]] + r[1] * (*prim_grad_y[iVal])[tempJ[0]]);
				riemann_para.primR(iVal) = (*prim[iVal])[iNode] + 0.5 * (*limiter_coef[iVal])[iNode] * (r[0] * (*prim_grad_x[iVal])[iNode] + r[1] * (*prim_grad_y[iVal])[iNode]);
			}
			riemannSolver_->Solver(riemann_para);
			for (int iVal = 0; iVal < GetNumberOfEquations(); ++iVal)
			{
				(*res[iVal])[iNode] -= riemann_para.flux[iVal] / jacobi;
				if (isnan((*res[iVal])[iNode]) || isinf((*res[iVal])[iNode]))
					ZaranLog::error("inode={},NAN in Residual!", iNode);
				(*cons[iVal])[iNode] -= (*res[iVal])[iNode] * dt[iNode] * (*coord_trans_coef[32])[iNode];
				(*res[iVal])[iNode] = 0;
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
		auto& prim = m_Primitive;
		auto& cons = m_Conservative;
		auto& prim_grad_x = m_PrimGradX;
		auto& prim_grad_y = m_PrimGradY;
		auto& res = m_Residual;
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
				(*prim[iVal])[iNode] = 0;
				(*res[iVal])[iNode] = 0;
				for (int iNeib = 0; iNeib < neiborNode.size(); ++iNeib)
				{
					(*prim[iVal])[iNode] += (*prim[iVal])[neiborNode[iNeib]] * weight[iNeib];
				}
			}
			Primitive2Conservative((*prim[0])[iNode], (*prim[1])[iNode], (*prim[2])[iNode], (*prim[3])[iNode], (*prim[4])[iNode], (*cons[0])[iNode], (*cons[1])[iNode], (*cons[2])[iNode], (*cons[3])[iNode], (*cons[4])[iNode]);
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
		auto& prim = m_Primitive;
		auto& cons = m_Conservative;
		auto& prim_grad_x = m_PrimGradX;
		auto& prim_grad_y = m_PrimGradY;
		auto& res = m_Residual;
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
				(*prim[iVal])[iNode] = prim_ideal(iVal);
				(*res[iVal])[iNode] = 0;
			}
			Primitive2Conservative((*prim[0])[iNode], (*prim[1])[iNode], (*prim[2])[iNode], (*prim[3])[iNode], (*prim[4])[iNode], (*cons[0])[iNode], (*cons[1])[iNode], (*cons[2])[iNode], (*cons[3])[iNode], (*cons[4])[iNode]);
		}
	}
	void Solver_NS_2D::ViscousFlux()
	{

	}

	void Solver_NS_2D::SourceFlux()
	{
	}

}