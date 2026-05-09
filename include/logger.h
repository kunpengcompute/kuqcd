/*
 * Copyright (c) 2026 Huawei Technologies Co., Ltd. All Rights Reserved.
 *
 * KuQCD is licensed under Mulan PSL v2.
 * You can use this software according to the terms and conditions of the Mulan PSL v2.
 * You may obtain a copy of Mulan PSL v2 at:
 *        http://license.coscl.org.cn/MulanPSL2
 * THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
 * See the Mulan PSL v2 for more details.
 */
#ifndef LOGGER_H
#define LOGGER_H

#include <cassert>
#include <cstdio>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <ostream>
#include <sstream>
#include <stack>
#include <string>
#include "kuqcd_complex.h"
#include "spinor.h"
#include "su3.h"

namespace kuqcd {

// 日志级别枚举
enum class LogLevel {
    ALL,
    ALLOC,
    TRACE,
    DEBUG,
    INFO,
    WARN,
    ERROR,
    FATAL,
    OFF
};
static const char *LogLevelStr[] = {"ALL", "ALLOC", "TRACE", "DEBUG", "INFO", "WARN", "ERROR", "FATAL", "OFF"};

// 颜色定义命名空间
namespace Colors {
const std::string red("\033[0;31m");
const std::string redBold("\033[1;31m");
const std::string green("\033[0;32m");
const std::string greenBold("\033[1;32m");
const std::string yellow("\033[0;33m");
const std::string yellowBold("\033[1;33m");
const std::string cyan("\033[0;36m");
const std::string cyanBold("\033[1;36m");
const std::string magenta("\033[0;35m");
const std::string magentaBold("\033[1;35m");
const std::string blue("\033[0;34m");
const std::string blueBold("\033[1;34m");
const std::string reset("\033[0m");
} // namespace Colors

// 字符串连接辅助函数
template <typename... Args>
inline std::string sjoin(Args &&...args) noexcept
{
    std::ostringstream msg;
    (msg << ... << args);
    return msg.str();
}

// 时间戳函数
inline std::string timeStamp()
{
    time_t rawtime;
    struct tm *timeinfo = nullptr;
    char buffer[80];

    time(&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer, sizeof(buffer), "%Y-%m-%d %H:%M:%S ", timeinfo);
    return std::string(buffer);
}

// Logger类
class Logger {
private:
    LogLevel log_level;
    bool colorized_output;
    std::stack<LogLevel> verbosity_stack;
    std::ostream &out;
    std::string addPrefix;

public:
    Logger(LogLevel loglevel = LogLevel::ALL, bool colorized_output = true, std::ostream &out = std::cout,
           std::string addPrefix = "");

    void SetVerbosity(LogLevel v);
    LogLevel GetVerbosity();

    void PushVerbosity(LogLevel verbosity);
    void PopVerbosity();
    void SetColoredOutput(bool val);
    void SetAdditionalPrefix(std::string add);
    std::string GetAdditionalPrefix();

    // 核心消息输出模板函数
    template <LogLevel level, typename... Args>
    std::string Message(Args &&...args)
    {
        std::ostringstream prefix, loginfo, postfix;

        bool resetColor = false;
        if (colorized_output && level == LogLevel::WARN) {
            prefix << Colors::yellow;
            resetColor = true;
        }
        if (colorized_output && level == LogLevel::ERROR) {
            prefix << Colors::red;
            resetColor = true;
        }
        if (colorized_output && level == LogLevel::FATAL) {
            prefix << Colors::redBold;
            resetColor = true;
        }

        loginfo << "# " << timeStamp() << LogLevelStr[static_cast<int>(level)] << ": ";

        std::string msg = sjoin(std::forward<Args>(args)...);

        if (resetColor)
            postfix << Colors::reset;

        postfix << std::resetiosflags(std::ios_base::floatfield | std::ios_base::basefield |
                                      std::ios_base::adjustfield | std::ios_base::uppercase | std::ios_base::showpos |
                                      std::ios_base::showpoint | std::ios_base::showbase | std::ios_base::boolalpha)
                << std::endl;

        if (level >= log_level) {
            if (addPrefix == "[ROOT] ") {
                if (GetRank() == 0) {
                    out << prefix.str() << loginfo.str() << addPrefix << msg << postfix.str();
                }
            } else {
                out << prefix.str() << loginfo.str() << "[" << std::setw(4) << std::setfill(' ') << GetRank() << "]"
                    << msg << postfix.str();
            }
        }

        return prefix.str() + addPrefix + msg + postfix.str();
    }

    // 各级别日志输出函数
    template <typename... Args>
    std::string Info(Args &&...args)
    {
        return Message<LogLevel::INFO>(std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::string Trace(Args &&...args)
    {
        return Message<LogLevel::TRACE>(std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::string Alloc(Args &&...args)
    {
        return Message<LogLevel::ALLOC>(std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::string Debug(Args &&...args)
    {
        return Message<LogLevel::DEBUG>(std::forward<Args>(args)...);
    }

    template <typename... Args>
    std::string Warn(Args &&...args)
    {
        return Message<LogLevel::WARN>(std::forward<Args>(args)...);
    }

    /*! 当出现错误但程序可以继续运行时使用
     * 例如：测试给出错误结果
     */
    template <typename... Args>
    std::string Error(Args &&...args)
    {
        return Message<LogLevel::ERROR>(std::forward<Args>(args)...);
    }

    /*! 当程序应该终止时使用（即抛出运行时异常）
     * 例如：存在冲突的输入参数
     */
    template <typename... Args>
    std::string Fatal(Args &&...args)
    {
        return Message<LogLevel::FATAL>(std::forward<Args>(args)...);
    }
};

// 全局logger实例声明
extern Logger stdLogger;
extern Logger rootLogger;

} // namespace kuqcd

// 为kuqcd类型提供流输出操作符重载
namespace kuqcd {

// Complex类的输出操作符
template <typename T, int Nbatch>
std::ostream &operator<<(std::ostream &os, const Complex<T, Nbatch> &c)
{
    if (Nbatch == 1) {
        if (c.Imag() >= 0) {
            os << c.Real() << "+" << c.Imag() << "i";
        } else {
            os << c.Real() << c.Imag() << "i";
        }
    } else {
        os << "Complex[" << Nbatch << "](";
        for (int i = 0; i < Nbatch; ++i) {
            if (i > 0)
                os << ", ";
            if (c.Imag(i) >= 0) {
                os << c.Real(i) << "+" << c.Imag(i) << "i";
            } else {
                os << c.Real(i) << c.Imag(i) << "i";
            }
        }
        os << ")";
    }
    return os;
}

// ColorSpinor类的输出操作符
template <typename T, int Nbatch>
std::ostream &operator<<(std::ostream &os, const ColorSpinor<T, Nbatch> &cs)
{
    os << "ColorSpinor(";
    os << "red=" << cs.red() << ", ";
    os << "green=" << cs.green() << ", ";
    os << "blue=" << cs.blue();
    os << ")";
    return os;
}

// SU3类的输出操作符
template <typename T, int Nbatch, KuQCDReconstructType Recon>
std::ostream &operator<<(std::ostream &os, const SU3<T, Nbatch, Recon> &su3)
{
    os << "SU3[" << static_cast<int>(Recon) << "](\n";

    // 根据压缩方案输出不同的信息
    if constexpr (Recon == KUQCD_GAUGE_RECON_NO) {
        // 完整3x3矩阵
        for (int i = 0; i < 3; ++i) {
            os << "  [";
            for (int j = 0; j < 3; ++j) {
                if (j > 0)
                    os << ", ";
                os << su3(i, j);
            }
            os << "]";
            if (i < 2)
                os << "\n";
        }
    } else {
        // 压缩格式，只显示存储的元素
        os << "  [compressed format]";
    }

    os << "\n)";
    return os;
}

} // namespace kuqcd
#endif
