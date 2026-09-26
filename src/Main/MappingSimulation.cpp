#include "MappingSimulation.h"
#include "GlobalData.h"
#include "Log.h"
#include "ZaranError.h"
#include "Visual.h"
#include "File.h"
#include "DEMFieldData.h"
#include "ReadDEMParticle.h"

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace zaran
{

namespace
{
	double ReadDoubleOr(const char* key, double fallback)
	{
		return GlobalData::IsExist(key) ? GlobalData::GetDouble(key) : fallback;
	}

	int ReadIntOr(const char* key, int fallback)
	{
		return GlobalData::IsExist(key) ? GlobalData::GetInt(key) : fallback;
	}

	bool ReadBoolOr(const char* key, bool fallback)
	{
		return GlobalData::IsExist(key) ? GlobalData::GetBool(key) : fallback;
	}

	std::string TrimNumber(double v, int digits = 3)
	{
		char buf[64];
		std::snprintf(buf, sizeof(buf), "%.*g", digits, v);
		return std::string(buf);
	}

	/// @brief 整块落盘（文本模式，与其它 ASCII 输出一致）
	bool WriteText(const std::string& path, const std::string& content)
	{
		std::ofstream out(path);
		if (!out.is_open())
		{
			Log::warn("MappingSimulation: cannot open '{}' for writing", path);
			return false;
		}
		out.write(content.data(), static_cast<std::streamsize>(content.size()));
		return out.good();
	}
}

MappingSimulation::MappingSimulation(std::string work_dir)
	: m_work_dir(std::move(work_dir))
{
}

void MappingSimulation::ReadParameters()
{
	m_result_folder = GlobalData::IsExist("output.result_folder")
		? GlobalData::GetString("output.result_folder") : "result";

	m_particle_file = GlobalData::IsExist("mapping.particle_file")
		? GlobalData::GetString("mapping.particle_file") : "particles.csv";

	// 计算域：给的是**物理域**（含两端半格），与 DEM 的 box 口径一致
	m_lo[0] = ReadDoubleOr("mapping.x_min", -0.03);
	m_hi[0] = ReadDoubleOr("mapping.x_max", 0.03);
	m_lo[1] = ReadDoubleOr("mapping.y_min", -0.03);
	m_hi[1] = ReadDoubleOr("mapping.y_max", 0.03);
	m_lo[2] = ReadDoubleOr("mapping.z_min", -0.03);
	m_hi[2] = ReadDoubleOr("mapping.z_max", 0.03);

	m_dim = ReadIntOr("mapping.dim", 3);
	if (m_dim != 2 && m_dim != 3)
	{
		throw ZaranError("MappingSimulation: [mapping] dim must be 2 or 3");
	}

	m_n0 = ReadIntOr("mapping.n0", 4);
	m_refine = ReadIntOr("mapping.refine", 2);
	m_levels = ReadIntOr("mapping.n_levels", 4);
	m_method = GlobalData::IsExist("mapping.method")
		? GlobalData::GetString("mapping.method") : "both";
	m_sub = ReadIntOr("mapping.n_sub", 2);
	m_smooth_passes = ReadIntOr("mapping.smooth_passes", 0);
	m_smooth_blend = ReadDoubleOr("mapping.smooth_blend", 1.0);
	m_write_particles = ReadBoolOr("mapping.write_particles", true);

	if (m_n0 < 1 || m_refine < 1 || m_levels < 1)
	{
		throw ZaranError("MappingSimulation: n0/refine/n_levels must be positive");
	}
	if (m_method != "subcell" && m_method != "exact" && m_method != "both")
	{
		throw ZaranError("MappingSimulation: [mapping] method must be subcell | exact | both");
	}
	for (int d = 0; d < m_dim; ++d)
	{
		if (!(m_hi[d] > m_lo[d]))
		{
			throw ZaranError("MappingSimulation: [mapping] box has a non-positive extent");
		}
	}

	Log::info("MappingSimulation: dim={}, box=[{},{}]x[{},{}]x[{},{}], levels={} (n0={}, refine={}), "
		"method={}, n_sub={}, smoothing={}x{:.2g}",
		m_dim, m_lo[0], m_hi[0], m_lo[1], m_hi[1], m_lo[2], m_hi[2],
		m_levels, m_n0, m_refine, m_method, m_sub, m_smooth_passes, m_smooth_blend);
}

void MappingSimulation::ReadParticles()
{
	const std::string path = m_work_dir + "/" + m_particle_file;
	if (!IsFileExist(path))
	{
		throw ZaranError("MappingSimulation: particle file not found: " + path);
	}
	// 复用 DEM 的读取器 ⇒ 输入格式与 DEM 算例完全一致
	ReadDEMParticle reader;
	std::vector<DEMParticle> raw;
	reader.ReadCSV(path, raw);
	if (raw.empty())
	{
		throw ZaranError("MappingSimulation: particle file has no particle: " + path);
	}

	m_particles.clear();
	m_particles.reserve(raw.size());
	m_particle_volume = 0.0;
	double min_radius = 1.0e30;
	m_outside_count = 0;

	for (const auto& p : raw)
	{
		if (!(p.radius > 0.0))
		{
			continue;
		}
		ParticleGridMapper::Particle q;
		q.x = p.pos.x();
		q.y = p.pos.y();
		q.z = p.pos.z();
		q.radius = p.radius;
		m_particles.push_back(q);

		const double vp = ParticleGridMapper::ParticleVolume(m_dim, q.radius);
		m_particle_volume += vp;
		min_radius = std::min(min_radius, q.radius);

		// 中心落在域外的粒子会在 BuildWeights 里被丢弃 ⇒ 体积账立刻不平。
		// 这里显式统计并在日志里点名，避免把"几何摆放错误"误读成"映射器 bug"。
		if (q.x < m_lo[0] || q.x > m_hi[0] || q.y < m_lo[1] || q.y > m_hi[1]
			|| (m_dim >= 3 && (q.z < m_lo[2] || q.z > m_hi[2])))
		{
			++m_outside_count;
		}
	}
	if (m_particles.empty())
	{
		throw ZaranError("MappingSimulation: no particle with a positive radius");
	}

	// 代表性粒径：用最小半径（映射最苛刻的那个），与文献判据 Δx ≥ 3–5 d_p 同口径
	m_dp = 2.0 * min_radius;

	Log::info("MappingSimulation: {} particles read from '{}'，ΣV_p={:.6e} m^3，d_p(min)={:.6e} m，"
		"域外粒子={}",
		m_particles.size(), m_particle_file, m_particle_volume, m_dp, m_outside_count);
	if (m_outside_count > 0)
	{
		Log::warn("MappingSimulation: {} 个粒子的中心在计算域外，这部分体积不会落到网格上，"
			"体积守恒残差会因此偏大", m_outside_count);
	}
}

ParticleGridMapper::GridSpec MappingSimulation::MakeGrid(int level) const
{
	ParticleGridMapper::GridSpec g;
	g.dim = m_dim;
	long long n = m_n0;
	for (int k = 0; k < level; ++k)
	{
		n *= m_refine;
	}
	g.ni = static_cast<int>(n);
	g.nj = (m_dim >= 3) ? static_cast<int>(n) : 1;
	g.nk = (m_dim >= 3) ? static_cast<int>(n) : 1;

	g.dx = (m_hi[0] - m_lo[0]) / g.ni;
	g.dy = (m_hi[1] - m_lo[1]) / g.nj;
	g.dz = (m_dim >= 3) ? (m_hi[2] - m_lo[2]) / g.nk : 1.0;

	// 节点型控制体：首节点 = 第 0 个单元的中心 ⇒ 两端各留半格
	g.x0 = m_lo[0] + 0.5 * g.dx;
	g.y0 = m_lo[1] + 0.5 * g.dy;
	g.z0 = (m_dim >= 3) ? (m_lo[2] + 0.5 * g.dz) : 0.0;
	return g;
}

void MappingSimulation::WriteFrame(int frame_index, int level, ParticleGridMapper::Method method,
	const ParticleGridMapper::GridSpec& grid,
	const dynamic_array<double>& alpha,
	const ParticleGridMapper::Stats& stats,
	bool append_header)
{
	Visual vis;
	const std::string dir = m_work_dir + "/" + m_result_folder;

	if (append_header && m_write_particles)
	{
		// 粒子点云：复用 DEM 的 Tecplot 写出器（点云与 α 云图可在 Tecplot 里叠加）
		DEMFieldData dem;
		for (size_t i = 0; i < m_particles.size(); ++i)
		{
			DEMParticle p;
			p.id = static_cast<index_type>(i);
			p.radius = m_particles[i].radius;
			p.pos = Eigen::Vector3d(m_particles[i].x, m_particles[i].y, m_particles[i].z);
			dem.AddParticle(p);
		}
		vis.WriteParticleTecplotASCII(dem, dir + "/particles.dat", std::string(), 0.0, false);
	}

	const std::string method_name = (method == ParticleGridMapper::Method::Exact) ? "exact" : "subcell";
	const double dx_over_dp = (m_dp > 0.0) ? grid.dx / m_dp : 0.0;

	// 帧名带上方法/网格/Δx·d_p⁻¹：Tecplot 的 Zone 列表里一眼就能分辨
	std::string zone = "alpha " + method_name
		+ " n=" + std::to_string(grid.ni)
		+ " dx/dp=" + TrimNumber(dx_over_dp);
	if (m_smooth_passes > 0)
	{
		zone += " smooth" + std::to_string(m_smooth_passes);
	}

	const std::string file_name = "alpha_" + std::to_string(frame_index) + ".dat";
	std::vector<Visual::UniformScalarField> fields;
	Visual::UniformScalarField f;
	f.name = "Alpha";
	f.values = &alpha;
	fields.push_back(f);

	vis.WriteUniformGridScalarTecplotASCII(dir + "/" + file_name,
		"Zaran3 solid volume fraction (particle -> grid)", zone,
		grid.ni, grid.nj, grid.nk, grid.x0, grid.y0, grid.z0,
		grid.dx, grid.dy, grid.dz, fields, static_cast<double>(level));

	FrameRecord rec;
	rec.index = frame_index;
	rec.method = method_name;
	rec.n_cell = grid.ni;
	rec.dx = grid.dx;
	rec.dx_over_dp = dx_over_dp;
	rec.residual = (m_particle_volume > 0.0)
		? std::fabs(stats.assigned_volume - m_particle_volume) / m_particle_volume
		: 0.0;
	rec.alpha_min = stats.alpha_min;
	rec.alpha_max = stats.alpha_max;
	rec.max_jump = stats.max_jump;
	rec.smooth_shift = stats.smoothing_shift;
	rec.file = file_name;
	m_records.push_back(rec);
}

void MappingSimulation::WriteReport() const
{
	const std::string path = m_work_dir + "/" + m_result_folder + "/mapping_report.csv";
	std::string buf;
	buf.reserve(m_records.size() * 128 + 256);
	buf += "frame,method,n_cell,dx,dx_over_dp,residual,alpha_min,alpha_max,"
		"max_jump,smooth_shift,l1_vs_exact,file\n";
	char line[320];
	for (const auto& r : m_records)
	{
		std::snprintf(line, sizeof(line),
			"%d,%s,%d,%.12e,%.6f,%.12e,%.12e,%.12e,%.12e,%.12e,%.12e,%s\n",
			r.index, r.method.c_str(), r.n_cell, r.dx, r.dx_over_dp, r.residual,
			r.alpha_min, r.alpha_max, r.max_jump, r.smooth_shift, r.l1_vs_exact,
			r.file.c_str());
		buf += line;
	}
	WriteText(path, buf);
	Log::info("MappingSimulation: report written to '{}'", path);
}

void MappingSimulation::Run()
{
	ReadParameters();
	ReadParticles();

	if (m_write_particles)
	{
		Log::info("MappingSimulation: 粒子点云会在首帧写出为 '{}/particles.dat'", m_result_folder);
	}

	const bool do_subcell = (m_method == "subcell" || m_method == "both");
	const bool do_exact = (m_method == "exact" || m_method == "both");

	Log::info("");
	Log::info("{:>5} {:>8} {:>11} {:>14} {:>14} {:>10} {:>10}",
		"frame", "method", "n_cell", "dx/dp", "ΣαV/Vp-1", "a_max", "max_jump");

	int frame_index = 0;
	for (int level = 0; level < m_levels; ++level)
	{
		const ParticleGridMapper::GridSpec grid = MakeGrid(level);

		dynamic_array<double> alpha_exact;
		ParticleGridMapper::Stats stats_exact;
		if (do_exact)
		{
			ParticleGridMapper mapper;
			mapper.Init(grid);
			mapper.SetMethod(ParticleGridMapper::Method::Exact);
			mapper.SetSmoothing(m_smooth_passes, m_smooth_blend);
			mapper.ComputeSolidFraction(m_particles, alpha_exact);
			stats_exact = mapper.GetStats();   // 按值留一份；下面还要用来对比

			WriteFrame(frame_index, level, ParticleGridMapper::Method::Exact, grid,
				alpha_exact, stats_exact, frame_index == 0);
			const FrameRecord& r = m_records.back();
			Log::info("{:>5} {:>8} {:>11} {:>14.4f} {:>14.3e} {:>10.4f} {:>10.4f}",
				r.index, r.method, r.n_cell, r.dx_over_dp, r.residual, r.alpha_max, r.max_jump);
			++frame_index;
		}

		if (do_subcell)
		{
			ParticleGridMapper mapper;
			mapper.Init(grid);
			mapper.SetMethod(ParticleGridMapper::Method::SubCell);
			mapper.SetSubdivisions(m_sub);
			mapper.SetSmoothing(m_smooth_passes, m_smooth_blend);
			dynamic_array<double> alpha;
			mapper.ComputeSolidFraction(m_particles, alpha);
			const ParticleGridMapper::Stats stats = mapper.GetStats();

			WriteFrame(frame_index, level, ParticleGridMapper::Method::SubCell, grid,
				alpha, stats, frame_index == 0);

			// 与同级 Exact 对比：这才是"SubCell 够不够用"的直接证据
			double l1 = -1.0;
			if (do_exact && alpha_exact.size() == alpha.size() && !alpha.empty())
			{
				double sum = 0.0;
				for (size_t c = 0; c < alpha.size(); ++c)
				{
					sum += std::fabs(alpha[c] - alpha_exact[c]);
				}
				l1 = sum / static_cast<double>(alpha.size());
			}
			m_records.back().l1_vs_exact = l1;

			const FrameRecord& r = m_records.back();
			Log::info("{:>5} {:>8} {:>11} {:>14.4f} {:>14.3e} {:>10.4f} {:>10.4f}",
				r.index, r.method, r.n_cell, r.dx_over_dp, r.residual, r.alpha_max, r.max_jump);
			++frame_index;
		}
	}

	WriteReport();

	// ---- 汇总：把 Δx/d_p 的方法选择准则当场量出来 ----
	Log::info("");
	Log::info("Δx/d_p 判据（L1 = 与同级 Exact 结果的单元平均绝对偏差）：");
	Log::info("{:>11} {:>14} {:>12} {:>12}", "dx/dp", "L1(sub-exact)", "a_max(sub)", "a_max(exact)");
	for (const auto& sub : m_records)
	{
		if (sub.method != "subcell" || sub.l1_vs_exact < 0.0)
		{
			continue;
		}
		double exact_max = 0.0;
		for (const auto& e : m_records)
		{
			if (e.method == "exact" && e.n_cell == sub.n_cell)
			{
				exact_max = e.alpha_max;
			}
		}
		Log::info("{:>11.4f} {:>14.3e} {:>12.4f} {:>12.4f}",
			sub.dx_over_dp, sub.l1_vs_exact, sub.alpha_max, exact_max);
	}
	Log::info("判据：Δx/d_p ≳ 3 时 SubCell(n_sub={}) 已足够；1–3 需要加大 n_sub 或改用 exact；",
		m_sub);
	Log::info("      < 1 时 SubCell 会把整块体积压进少数单元（α 出现尖峰），必须用 exact；");
	Log::info("      < 0.5 时体积分数法本身失效，需要转点源法。");

	const std::string dir = m_work_dir + "/" + m_result_folder;
	Log::info("MappingSimulation finished. {} frames written to '{}'", m_records.size(), dir);
}

} // namespace zaran
