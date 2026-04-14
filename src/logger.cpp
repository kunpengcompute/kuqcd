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
#include "logger.h"

namespace kuqcd {

// Logger类构造函数实现
Logger::Logger(LogLevel loglevel, bool colorized_output, std::ostream &out, std::string addPrefix)
    : log_level(loglevel), colorized_output(colorized_output), out(out), addPrefix(addPrefix)
{}

// 成员函数实现
void Logger::SetVerbosity(LogLevel v)
{
    log_level = v;
}

LogLevel Logger::GetVerbosity()
{
    return log_level;
}

void Logger::PushVerbosity(LogLevel verbosity)
{
    verbosity_stack.push(GetVerbosity());
    SetVerbosity(verbosity);
}

void Logger::PopVerbosity()
{
    SetVerbosity(verbosity_stack.top());
    verbosity_stack.pop();
}

void Logger::SetColoredOutput(bool val)
{
    colorized_output = val;
}

void Logger::SetAdditionalPrefix(std::string add)
{
    addPrefix = add;
}

std::string Logger::GetAdditionalPrefix()
{
    return addPrefix;
}

// 全局logger实例定义
Logger stdLogger(LogLevel::INFO, true, std::cout, "");
Logger rootLogger(LogLevel::INFO, true, std::cout, "[ROOT] ");

}  // namespace kuqcd