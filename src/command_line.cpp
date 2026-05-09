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
#include "command_line.h"
#include <iostream>
#include <string>
#include <map>
#include <stdexcept>
#include <mpi.h>
#include "logger.h"

namespace kuqcd {

std::string Trim(const std::string &s)
{
    size_t first = s.find_first_not_of(" \t");
    if (first == std::string::npos)
        return "";
    size_t last = s.find_last_not_of(" \t");
    return s.substr(first, (last - first + 1));
}

void CmdLineParam::ReadParam()
{
    int flag;
    MPI_Initialized(&flag);
    if (!flag) {
        rootLogger.Error("MPI is not initialized!");
        return;
    }
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        std::map<std::string, int *> variables = {
            {"gridx", &gridx}, {"gridy", &gridy}, {"gridz", &gridz}, {"gridt", &gridt}};

        std::string line;
        while (getline(std::cin, line)) {
            size_t eq_pos = line.find('=');
            if (eq_pos == std::string::npos) {
                std::cerr << "ignore invalid lines (no =): " << line << std::endl;
                continue;
            }

            std::string key = Trim(line.substr(0, eq_pos));
            std::string value_str = Trim(line.substr(eq_pos + 1));

            auto it = variables.find(key);
            if (it == variables.end()) {
                std::cerr << "unknown variable" << key << "'" << std::endl;
                continue;
            }

            try {
                int value = stoi(value_str);
                *(it->second) = value;
            } catch (const std::invalid_argument &e) {
                std::cerr << "invalid number '" << value_str << "' corresponding variable '" << key << "'" << std::endl;
            } catch (const std::out_of_range &e) {
                std::cerr << "number out of range '" << value_str << "' corresponding variable '" << key << "'"
                          << std::endl;
            }
        }

        std::cout << "gridx = " << gridx << "  ";
        std::cout << "gridy = " << gridy << "  ";
        std::cout << "gridz = " << gridz << "  ";
        std::cout << "gridt = " << gridt << std::endl;
    }

    MPI_Bcast(&gridx, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gridy, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gridz, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&gridt, 1, MPI_INT, 0, MPI_COMM_WORLD);
}

} // namespace kuqcd
