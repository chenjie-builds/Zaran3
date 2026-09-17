/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file mapping_main.cpp
 * \brief 粒子→网格映射（ParticleGridMapper）的基准/验证驱动。
 * \author Chen Jie.
 *
 *  这是**纯几何工具**的验证程序，不参与主程序构建（独立 CMake target MappingBench）。
 *  输出全部是 `key=value` 行或 CSV，交给 benchmarks/two_phase_mapping/verify_mapping.py
 *  做判定。随机粒子云用自定义 LCG 生成，Python 侧用同样的 64 位整数运算复刻，
 *  保证两侧看到完全相同的粒子构型。
 *
 *  模式
 *  ----
 *   field     转储固相体积分数场 + 粒子列表（供 Python 与子采样参考解对比）
 *   summary   汇总各方法/细分数下的配分残差与统计量
 *   force     力散射与对偶插值：Σ_cell = Σ_particle、常数场插值回原值
 *   smooth    保守光顺：Σα 变化量与最大跳变的下降
 */

#include "ParticleGridMapper.h"
#include "Log.h"
#include "ZaranError.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>

using namespace zaran;

namespace
{
	// ----------------------------------------------------------------
	// 确定性 PRNG（splitmix 风格的 64 位 LCG；Python 侧逐位复刻）
	// ----------------------------------------------------------------
	inline unsigned long long NextState(unsigned long long s)
	{
		return s * 6364136223846793005ULL + 1442695040888963407ULL;
	}
	inline double NextUnit(unsigned long long& s)
	{
		s = NextState(s);
		return static_cast<double>(s >> 11) * (1.0 / 9007199254740992.0);
	}

	struct Options
	{
		std::string mode = "summary";
		std::string out;
		int ni = 40, nj = 40, nk = 40;
		int dim = 3;
		double lx = 1.0, ly = 1.0, lz = 1.0;
		double radius = 0.02;
		int np = 200;
		int sub = 2;
		int smooth = 0;
		double blend = 1.0;
		unsigned long long seed = 20260917ULL;
		std::string method = "subcell";
		/// @brief 给定则只放一个粒子在 (px,py,pz)（用于对称性/单个单元检验）
		double px = 0.0, py = 0.0, pz = 0.0;
		int single = 0;
	};

	Options Parse(int argc, char** argv)
	{
		Options o;
		for (int i = 1; i < argc; ++i)
		{
			const std::string a = argv[i];
			auto val = [&](const char* name) -> std::string
			{
				if (i + 1 >= argc)
				{
					throw ZaranError(std::string("missing value for ") + name);
				}
				return argv[++i];
			};
			if (a == "--mode") o.mode = val("--mode");
			else if (a == "--out") o.out = val("--out");
			else if (a == "--ni") o.ni = std::atoi(val("--ni").c_str());
			else if (a == "--nj") o.nj = std::atoi(val("--nj").c_str());
			else if (a == "--nk") o.nk = std::atoi(val("--nk").c_str());
			else if (a == "--dim") o.dim = std::atoi(val("--dim").c_str());
			else if (a == "--lx") o.lx = std::atof(val("--lx").c_str());
			else if (a == "--ly") o.ly = std::atof(val("--ly").c_str());
			else if (a == "--lz") o.lz = std::atof(val("--lz").c_str());
			else if (a == "--radius") o.radius = std::atof(val("--radius").c_str());
			else if (a == "--np") o.np = std::atoi(val("--np").c_str());
			else if (a == "--sub") o.sub = std::atoi(val("--sub").c_str());
			else if (a == "--smooth") o.smooth = std::atoi(val("--smooth").c_str());
			else if (a == "--blend") o.blend = std::atof(val("--blend").c_str());
			else if (a == "--seed") o.seed = std::strtoull(val("--seed").c_str(), nullptr, 10);
			else if (a == "--method") o.method = val("--method");
			else if (a == "--px") { o.px = std::atof(val("--px").c_str()); o.single = 1; }
			else if (a == "--py") { o.py = std::atof(val("--py").c_str()); o.single = 1; }
			else if (a == "--pz") { o.pz = std::atof(val("--pz").c_str()); o.single = 1; }
			else
			{
				throw ZaranError("unknown argument: " + a);
			}
		}
		return o;
	}

	ParticleGridMapper::GridSpec MakeGrid(const Options& o)
	{
		ParticleGridMapper::GridSpec g;
		g.ni = o.ni;
		g.nj = o.nj;
		g.nk = (o.dim >= 3) ? o.nk : 1;
		g.dim = o.dim;
		g.dx = o.lx / static_cast<double>(g.ni - 1);
		g.dy = o.ly / static_cast<double>(g.nj - 1);
		g.dz = (o.dim >= 3) ? (o.lz / static_cast<double>(g.nk - 1)) : 1.0;
		g.x0 = 0.0;
		g.y0 = 0.0;
		g.z0 = 0.0;
		return g;
	}

	/// @brief 粒子云：默认均匀随机；给定 --px/--py/--pz 时只放一个指定位置的粒子
	dynamic_array<ParticleGridMapper::Particle> MakeCloud(
		const ParticleGridMapper::GridSpec& g, const Options& o)
	{
		dynamic_array<ParticleGridMapper::Particle> ps;
		if (o.single)
		{
			ParticleGridMapper::Particle p;
			p.radius = o.radius;
			p.x = o.px;
			p.y = o.py;
			p.z = o.pz;
			ps.push_back(p);
			return ps;
		}
		ps.reserve(static_cast<size_t>(o.np));
		unsigned long long s = o.seed;
		const double lo_x = g.x0 - 0.5 * g.dx + o.radius;
		const double hi_x = g.x0 + (g.ni - 0.5) * g.dx - o.radius;
		const double lo_y = g.y0 - 0.5 * g.dy + o.radius;
		const double hi_y = g.y0 + (g.nj - 0.5) * g.dy - o.radius;
		const double lo_z = g.z0 - 0.5 * g.dz + o.radius;
		const double hi_z = g.z0 + (g.nk - 0.5) * g.dz - o.radius;
		for (int i = 0; i < o.np; ++i)
		{
			ParticleGridMapper::Particle p;
			p.radius = o.radius;
			p.x = lo_x + (hi_x - lo_x) * NextUnit(s);
			p.y = lo_y + (hi_y - lo_y) * NextUnit(s);
			p.z = (g.dim >= 3) ? (lo_z + (hi_z - lo_z) * NextUnit(s)) : 0.0;
			ps.push_back(p);
		}
		return ps;
	}

	void SetMethod(ParticleGridMapper& m, const Options& o)
	{
		m.SetMethod(o.method == "exact"
			? ParticleGridMapper::Method::Exact
			: ParticleGridMapper::Method::SubCell);
		m.SetSubdivisions(o.sub);
		m.SetSmoothing(o.smooth, o.blend);
	}

	// ----------------------------------------------------------------
	// summary：配分残差与统计
	// ----------------------------------------------------------------
	int RunSummary(const Options& o)
	{
		const auto g = MakeGrid(o);
		const auto ps = MakeCloud(g, o);
		const double cell_volume = g.CellVolume();
		double total_particle_volume = 0.0;
		for (const auto& p : ps)
		{
			total_particle_volume += (g.dim >= 3)
				? (4.0 / 3.0) * 3.14159265358979323846 * p.radius * p.radius * p.radius
				: 3.14159265358979323846 * p.radius * p.radius;
		}
		std::printf("grid=%d,%d,%d dim=%d dx=%.17g cell_volume=%.17g\n",
			g.ni, g.nj, g.nk, g.dim, g.dx, cell_volume);
		std::printf("particles=%d radius=%.17g total_particle_volume=%.17g\n",
			o.np, o.radius, total_particle_volume);

		struct Variant { const char* method; int sub; int smooth; };
		const Variant variants[] = {
			{ "subcell", 1, 0 }, { "subcell", 2, 0 }, { "subcell", 3, 0 },
			{ "subcell", 4, 0 }, { "subcell", 8, 0 }, { "exact", 0, 0 },
			{ "subcell", 2, 4 }, { "exact", 0, 4 },
		};
		for (const auto& v : variants)
		{
			Options vo = o;
			vo.method = v.method;
			vo.sub = v.sub;
			vo.smooth = v.smooth;
			vo.blend = 1.0;
			ParticleGridMapper mapper(g);
			SetMethod(mapper, vo);
			dynamic_array<double> alpha;
			const double assigned = mapper.ComputeSolidFraction(ps, alpha);
			const auto& st = mapper.GetStats();
			std::printf("variant method=%s sub=%d smooth=%d assigned_volume=%.17g "
				"rel_residual=%.17g alpha_min=%.17g alpha_max=%.17g max_jump=%.17g "
				"touched=%lld pairs=%lld smooth_shift=%.17g\n",
				v.method, v.sub, v.smooth, assigned,
				std::fabs(assigned - total_particle_volume) / total_particle_volume,
				st.alpha_min, st.alpha_max, st.max_jump,
				st.cells_touched, st.particle_cell_pairs, st.smoothing_shift);
		}
		return 0;
	}

	// ----------------------------------------------------------------
	// field：转储网格量 + 粒子列表
	// ----------------------------------------------------------------
	int RunField(const Options& o)
	{
		if (o.out.empty())
		{
			throw ZaranError("--mode field requires --out");
		}
		const auto g = MakeGrid(o);
		const auto ps = MakeCloud(g, o);
		ParticleGridMapper mapper(g);
		SetMethod(mapper, o);
		dynamic_array<double> alpha;
		const double assigned = mapper.ComputeSolidFraction(ps, alpha);
		const auto& st = mapper.GetStats();

		{
			std::ofstream f(o.out);
			f.precision(17);
			f << "i,j,k,x,y,z,alpha\n";
			for (int k = 0; k < g.nk; ++k)
			{
				for (int j = 0; j < g.nj; ++j)
				{
					for (int i = 0; i < g.ni; ++i)
					{
						const long long c = mapper.Idx(i, j, k);
						f << i << "," << j << "," << k << ","
							<< (g.x0 + i * g.dx) << "," << (g.y0 + j * g.dy) << ","
							<< (g.z0 + k * g.dz) << "," << alpha[static_cast<size_t>(c)] << "\n";
					}
				}
			}
		}
		{
			std::ofstream f(o.out + ".particles");
			f.precision(17);
			f << "x,y,z,r\n";
			for (const auto& p : ps)
			{
				f << p.x << "," << p.y << "," << p.z << "," << p.radius << "\n";
			}
		}
		std::printf("grid=%d,%d,%d dim=%d dx=%.17g dy=%.17g dz=%.17g cell_volume=%.17g\n",
			g.ni, g.nj, g.nk, g.dim, g.dx, g.dy, g.dz, g.CellVolume());
		std::printf("method=%s sub=%d smooth=%d alpha_min=%.17g alpha_max=%.17g "
			"max_jump=%.17g assigned=%.17g\n",
			o.method.c_str(), o.sub, o.smooth, st.alpha_min, st.alpha_max,
			st.max_jump, assigned);
		return 0;
	}

	// ----------------------------------------------------------------
	// force：力散射守恒 + 对偶插值
	// ----------------------------------------------------------------
	int RunForce(const Options& o)
	{
		const auto g = MakeGrid(o);
		const auto ps = MakeCloud(g, o);
		ParticleGridMapper mapper(g);
		SetMethod(mapper, o);

		// 每个粒子一个随机力
		dynamic_array<double> fp(ps.size() * 3, 0.0);
		unsigned long long s = o.seed ^ 0x9e3779b97f4a7c15ULL;
		for (size_t i = 0; i < fp.size(); ++i)
		{
			fp[i] = NextUnit(s) * 2.0 - 1.0;
		}
		dynamic_array<double> fc;
		mapper.ScatterParticleVector(ps, fp, fc);

		double sum_p[3] = { 0, 0, 0 }, sum_c[3] = { 0, 0, 0 };
		for (size_t i = 0; i < ps.size(); ++i)
		{
			for (int c = 0; c < 3; ++c) sum_p[c] += fp[i * 3 + c];
		}
		for (size_t i = 0; i < fc.size() / 3; ++i)
		{
			for (int c = 0; c < 3; ++c) sum_c[c] += fc[i * 3 + c];
		}
		double rel = 0.0;
		for (int c = 0; c < 3; ++c)
		{
			rel = std::max(rel, std::fabs(sum_c[c] - sum_p[c]) / std::max(std::fabs(sum_p[c]), 1e-300));
		}

		// 对偶性：常数场 gather 回常数
		dynamic_array<double> constant(fc.size(), 0.0);
		for (size_t i = 0; i < constant.size() / 3; ++i)
		{
			constant[i * 3 + 0] = 1.25;
			constant[i * 3 + 1] = -0.5;
			constant[i * 3 + 2] = 3.0;
		}
		dynamic_array<double> fg;
		mapper.GatherToParticle(ps, constant, fg);
		double dual = 0.0;
		for (size_t i = 0; i < ps.size(); ++i)
		{
			dual = std::max(dual, std::fabs(fg[i * 3 + 0] - 1.25));
			dual = std::max(dual, std::fabs(fg[i * 3 + 1] + 0.5));
			dual = std::max(dual, std::fabs(fg[i * 3 + 2] - 3.0));
		}
		std::printf("method=%s sub=%d particles=%d\n", o.method.c_str(), o.sub, o.np);
		std::printf("scatter_sum_particle=%.17g,%.17g,%.17g\n", sum_p[0], sum_p[1], sum_p[2]);
		std::printf("scatter_sum_cell=%.17g,%.17g,%.17g\n", sum_c[0], sum_c[1], sum_c[2]);
		std::printf("scatter_rel_residual=%.17g\n", rel);
		std::printf("gather_dual_max_error=%.17g\n", dual);
		return 0;
	}

	// ----------------------------------------------------------------
	// smooth：保守光顺
	// ----------------------------------------------------------------
	int RunSmooth(const Options& o)
	{
		const auto g = MakeGrid(o);
		const auto ps = MakeCloud(g, o);
		ParticleGridMapper plain(g);
		Options vo = o;
		vo.smooth = 0;
		SetMethod(plain, vo);
		dynamic_array<double> raw;
		const double assigned = plain.ComputeSolidFraction(ps, raw);
		double sum_raw = 0.0;
		for (double a : raw) sum_raw += a;

		std::printf("particles=%d radius=%.17g dx=%.17g dx_over_d=%.17g\n",
			o.np, o.radius, g.dx, g.dx / (2.0 * o.radius));
		std::printf("raw sum=%.17g max_jump=%.17g assigned=%.17g\n",
			sum_raw, plain.GetStats().max_jump, assigned);
		for (int passes : { 1, 2, 4, 8, 16 })
		{
			ParticleGridMapper sm(g);
			vo.smooth = passes;
			vo.blend = 1.0;
			SetMethod(sm, vo);
			dynamic_array<double> a;
			sm.ComputeSolidFraction(ps, a);
			double s1 = 0.0, amin = 1e300, amax = -1e300;
			for (double v : a)
			{
				s1 += v;
				amin = std::min(amin, v);
				amax = std::max(amax, v);
			}
			std::printf("smooth passes=%d sum=%.17g sum_rel_shift=%.17g max_jump=%.17g "
				"alpha_min=%.17g alpha_max=%.17g\n",
				passes, s1, std::fabs(s1 - sum_raw) / sum_raw, sm.GetStats().max_jump,
				amin, amax);
		}
		return 0;
	}
}

int main(int argc, char** argv)
{
	try
	{
		// 本驱动用 stdout 输出机器可读的 key=value，日志必须静音，否则会污染解析
		spdlog::set_level(spdlog::level::off);
		// stdout 无缓冲：即使中途崩溃也保证已经算出的结果能被脚本看到
		std::setvbuf(stdout, nullptr, _IONBF, 0);
		const Options o = Parse(argc, argv);
		if (o.mode == "summary") return RunSummary(o);
		if (o.mode == "field") return RunField(o);
		if (o.mode == "force") return RunForce(o);
		if (o.mode == "smooth") return RunSmooth(o);
		std::cerr << "unknown mode: " << o.mode << std::endl;
		return 2;
	}
	catch (const std::exception& e)
	{
		std::cerr << "MappingBench error: " << e.what() << std::endl;
		return 1;
	}
}
