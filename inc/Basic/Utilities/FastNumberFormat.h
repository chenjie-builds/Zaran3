/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file FastNumberFormat.h
 * \brief 轻量快速数值格式化：把整数/浮点追加到 std::string。
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 */
#pragma once
#include <charconv>
#include <string>
#include <type_traits>

namespace zaran
{
    /// @brief 追加整数到字符串缓冲。
    /// 相比 `ostream <<`，省去 locale/虚函数/流状态开销。
    template <typename T, typename = std::enable_if_t<std::is_integral_v<T>>>
    inline void AppendNumber(std::string& s, T v)
    {
        char buf[32];
        const auto r = std::to_chars(buf, buf + sizeof(buf), v);
        s.append(buf, static_cast<std::size_t>(r.ptr - buf));
    }

    /// @brief 追加浮点数到字符串缓冲（general 格式，6 位有效数字，
    /// 与 `std::ostream` 默认精度一致，故输出文本与旧实现保持可比）。
    inline void AppendNumber(std::string& s, double v)
    {
        char buf[40];
        const auto r = std::to_chars(buf, buf + sizeof(buf), v,
                                     std::chars_format::general, 6);
        s.append(buf, static_cast<std::size_t>(r.ptr - buf));
    }

    inline void AppendNumber(std::string& s, float v)
    {
        AppendNumber(s, static_cast<double>(v));
    }
} // namespace zaran
