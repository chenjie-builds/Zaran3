#include "GradWLSQ.h"
#include "Log.h"
using namespace zaran;
zaran::GradWLSQ::GradWLSQ(std::shared_ptr<GridFN>grid) : m_grid(grid)
{
	m_grid = grid;
	auto node = grid->GetNode();
	int node_num = node->GetCount();
	double r11, r12, r13, r22, r23, r23_a, r23_b, r33;
	double alpha1, alpha2, alpha3;
	double beta;
	double alpha[3];
	std::vector<double> dx, dy, dz, weight;
	m_omega = new double** [node_num];
	for (int iNode = 0; iNode < node_num; iNode++)
	{
		// if (node->GetType(iNode) != NodeType::inner)
		//     continue;
		r11 = r12 = r13 = r22 = r23 = r23_a = r23_b = r33 = 0.0;
		beta = 0.0;
		auto neighbors = node->GetNeighborNode(iNode);
		int neighbor_num = node->GetNeighborNodeNum(iNode);
		dx.resize(neighbor_num);
		dy.resize(neighbor_num);
		dz.resize(neighbor_num);
		weight.resize(neighbor_num);
		m_omega[iNode] = new double* [neighbor_num];
		for (int iNeigh = 0; iNeigh < neighbor_num; iNeigh++)
		{
			dx[iNeigh] = node->GetCoord(neighbors[iNeigh])[0] - node->GetCoord(iNode)[0];
			dy[iNeigh] = node->GetCoord(neighbors[iNeigh])[1] - node->GetCoord(iNode)[1];
			dz[iNeigh] = node->GetCoord(neighbors[iNeigh])[2] - node->GetCoord(iNode)[2];
			weight[iNeigh] = 1.0 / (dx[iNeigh] * dx[iNeigh] + dy[iNeigh] * dy[iNeigh] + dz[iNeigh] * dz[iNeigh]);
			weight[iNeigh] = sqrt(weight[iNeigh]);
			dx[iNeigh] *= weight[iNeigh];
			dy[iNeigh] *= weight[iNeigh];
			dz[iNeigh] *= weight[iNeigh];
			r11 += dx[iNeigh] * dx[iNeigh];
			r12 += dx[iNeigh] * dy[iNeigh];
			r22 += dy[iNeigh] * dy[iNeigh];
			r13 += dx[iNeigh] * dz[iNeigh];
			r23_a += dy[iNeigh] * dz[iNeigh];
			r23_b += dx[iNeigh] * dz[iNeigh];
			r33 += dz[iNeigh] * dz[iNeigh];
		}
		if (r11 >= 0.0)
			r11 = sqrt(r11);
		else
			r11 = 0.0;
		if (r11 != 0.0)
			r12 = r12 / r11;
		else
			r12 = 0.0;
		if (r22 - r12 * r12 >= 0.0)
			r22 = sqrt(r22 - r12 * r12);
		else
			r22 = 0.0;
		if (r11 != 0.0)
			r13 = r13 / r11;
		else
			r13 = 0.0;
		if (r22 != 0.0 && r11 != 0.0)
			r23 = r23_a / r22 - r23_b * r12 / (r11 * r22);
		else
			r23 = 0.0;
		if (r33 - r23 * r23 - r13 * r13 >= 0.0)
			r33 = sqrt(r33 - r23 * r23 - r13 * r13);
		else
			r33 = 0.0;
		if (r11 * r22 == 0.0)
			beta = 0.0;
		else
			beta = (r12 * r23 - r13 * r22) / (r11 * r22);
		for (int iNeigh = 0; iNeigh < neighbor_num; ++iNeigh)
		{
			m_omega[iNode][iNeigh] = new double[3];
			alpha[0] = dx[iNeigh] / (r11 * r11);
			alpha[1] = (dy[iNeigh] - r12 / r11 * dx[iNeigh]) / (r22 * r22);
			alpha[2] = (dz[iNeigh] - r23 / r22 * dy[iNeigh] + beta * dx[iNeigh]) / (r33 * r33);
			m_omega[iNode][iNeigh][0] = alpha[0] - r12 / r11 * alpha[1] + beta * alpha[2];
			m_omega[iNode][iNeigh][1] = alpha[1] - r23 / r22 * alpha[2];
			m_omega[iNode][iNeigh][2] = alpha[2];
			m_omega[iNode][iNeigh][0] *= weight[iNeigh];
			m_omega[iNode][iNeigh][1] *= weight[iNeigh];
			m_omega[iNode][iNeigh][2] *= weight[iNeigh];
			// if (isnan(m_omega[iNode][iNeigh][0]) || isnan(m_omega[iNode][iNeigh][1]) || isnan(m_omega[iNode][iNeigh][2]))
			// {
			//     Log::warn("GradWLSQ::GradWLSQ: omega is nan");
			// }
			if (std::isnan(m_omega[iNode][iNeigh][0]))
				m_omega[iNode][iNeigh][0] = 0.0;
			if (std::isnan(m_omega[iNode][iNeigh][1]))
				m_omega[iNode][iNeigh][1] = 0.0;
			if (std::isnan(m_omega[iNode][iNeigh][2]))
				m_omega[iNode][iNeigh][2] = 0.0;
		}
	}
}
GradWLSQ::~GradWLSQ()
{
	auto node = m_grid->GetNode();
	int node_num = node->GetCount();
	for (int iNode = 0; iNode < node_num; iNode++)
	{
		auto neighbors = node->GetNeighborNode(iNode);
		int neighbor_num = node->GetNeighborNodeNum(iNode);
		for (int iNeigh = 0; iNeigh < neighbor_num; iNeigh++)
		{
			delete[] m_omega[iNode][iNeigh];
		}
		delete[] m_omega[iNode];
	}
	delete[] m_omega;
}

void GradWLSQ::CalcGradient(std::shared_ptr<GridFN>grid, const double* data, double* grad_x, double* grad_y, double* grad_z)
{
	if (grid != m_grid)
	{
		Log::error("Gradient::CalcGradient: grid is not matched");
		return;
	}
	auto node = grid->GetNode();
	int node_num = node->GetCount();
#pragma omp parallel for
	for (int iNode = 0; iNode < node_num; iNode++)
	{
		if (node->GetType(iNode) != NodeType::inner)
			continue;
		auto neighbors = node->GetNeighborNode(iNode);
		int neighbor_num = node->GetNeighborNodeNum(iNode);
		grad_x[iNode] = grad_y[iNode] = grad_z[iNode] = 0.0;
		for (int iNeigh = 0; iNeigh < neighbor_num; iNeigh++)
		{
			int iNeighNode = neighbors[iNeigh];
			grad_x[iNode] += m_omega[iNode][iNeigh][0] * (data[iNeighNode] - data[iNode]);
			grad_y[iNode] += m_omega[iNode][iNeigh][1] * (data[iNeighNode] - data[iNode]);
			grad_z[iNode] += m_omega[iNode][iNeigh][2] * (data[iNeighNode] - data[iNode]);
		}
		if (std::isnan(grad_x[iNode]) || std::isnan(grad_y[iNode]) || std::isnan(grad_z[iNode]))
		{
			// Log::warn("GradWLSQ::CalcGradient: grad is nan");
			grad_x[iNode] = grad_y[iNode] = grad_z[iNode] = 0.0;
		}
	}
}