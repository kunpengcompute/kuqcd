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
#include "kuqcd_utils.h"

#ifdef USE_KUPL
#include <unistd.h>
#include "kupl.h"

int oob_allgather_callback(const void *sendbuf, void *recvbuf, int size, void *group, kupl_shm_datatype_t datatype)
{
    switch (datatype) {
        case KUPL_SHM_DATATYPE_CHAR:
            return MPI_Allgather(sendbuf, size, MPI_CHAR, recvbuf, size, MPI_CHAR, (MPI_Comm)group);
        case KUPL_SHM_DATATYPE_INT:
            return MPI_Allgather(sendbuf, size, MPI_INT, recvbuf, size, MPI_INT, (MPI_Comm)group);
        case KUPL_SHM_DATATYPE_LONG:
            return MPI_Allgather(sendbuf, size, MPI_LONG, recvbuf, size, MPI_LONG, (MPI_Comm)group);
        case KUPL_SHM_DATATYPE_FLOAT:
            return MPI_Allgather(sendbuf, size, MPI_FLOAT, recvbuf, size, MPI_FLOAT, (MPI_Comm)group);
        case KUPL_SHM_DATATYPE_DOUBLE:
            return MPI_Allgather(sendbuf, size, MPI_DOUBLE, recvbuf, size, MPI_DOUBLE, (MPI_Comm)group);
        default:
            fprintf(stderr, "not support datatype");
            return KUPL_ERROR;
    }
}
int oob_barrier_callback(void *group)
{
    return MPI_Barrier((MPI_Comm)group);
}
#endif

namespace kuqcd {
bool pathExists(const std::string &path)
{
    struct stat info;
    return stat(path.c_str(), &info) == 0 && S_ISDIR(info.st_mode);
}
bool numaHasHugePages(int nodeId)
{
    // 构建NUMA节点路径
    std::string basePath = "/sys/devices/system/node/node" + std::to_string(nodeId) + "/hugepages";

    // 检查节点目录是否存在
    if (!pathExists(basePath)) {
        printf("NUMA node %d is not exist\n", nodeId);
        return false;
    }

    DIR *dir = opendir(basePath.c_str());
    if (!dir) {
        perror("Failed to open directory");
        return false;
    }

    bool foundHugePages = false;
    struct dirent *entry;

    // 遍历所有大页尺寸目录
    while ((entry = readdir(dir)) != nullptr) {
        if (entry->d_type == DT_DIR && strncmp(entry->d_name, "hugepages-", 10) == 0) {
            std::string sizeDir = basePath + "/" + entry->d_name;
            std::string nrFile = sizeDir + "/nr_hugepages";

            // 读取nr_hugepages文件内容
            std::ifstream file(nrFile);
            if (file.is_open()) {
                int nrHugepages;
                file >> nrHugepages;
                if (nrHugepages > 0) {
                    foundHugePages = true;
                    break;
                }
            }
        }
    }
    closedir(dir);
    return foundHugePages;
}
} // namespace kuqcd
