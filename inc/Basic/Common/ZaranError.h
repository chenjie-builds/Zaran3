/**
 * Zaran	-	A Totally Automatic CFD Software
 * \file ZaranError.h
 * \brief Program-level fatal error exception, used to replace exit()/system("pause").
 * \author Chen Jie.
 *
 * \copyright Copyright (C) Since 2020, Chen Jie.
 * This file is part of Zaran.
 * All rights reserved. This software is proprietary and confidential.
 * Unauthorized copying, distribution, or use is strictly prohibited.
 */
#pragma once
#include <stdexcept>
#include <string>

namespace zaran
{
    /// @brief 程序级致命错误异常
    /// @details 用于取代 exit()/system("pause")，由 main() 顶层捕获后打印并以非零码退出
    class ZaranError : public std::runtime_error
    {
    public:
        explicit ZaranError(const std::string& msg) : std::runtime_error(msg) {}
    };
} // namespace zaran
