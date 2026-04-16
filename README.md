# KUQCD

## 🔥RELEASE-NOTES
- [2026/03] KuQCD 项目首次上线，支持HISQ dslash的Multishift CG求解器。

## 🚀概述

鲲鹏量子格点色动力学加速库（Kunpeng Unified Lattice Quantum Chromodynamics Library，以下简称KuQCD）是基于鲲鹏平台的格点量子色动力学加速库，采用C/C++实现，目前仅支持在鲲鹏920新型号上使用，依赖KUPL，支持毕昇编译器和gcc编译器。本加速库提供inverter等LQCD领域基础功能模块，可选不同求解算法和dslash，充分发挥鲲鹏处理器硬件特性，提供高性能基础接口。

## 📝版本配套

- 运行平台
    - 鲲鹏 920 专业版
- 系统规格
    - openEuler 22.03（LTS-SP4）AArch64

## ⚡️编译安装

### 1. [获取 HPCKit 软件包](https://www.hikunpeng.com/document/detail/zh/kunpenghpcs/hpckit/instg/KunpengHPCKit_install_007.html)
> https://www.hikunpeng.com/developer/hpc/hpckit-download
### 2. [安装 HPCKit](https://www.hikunpeng.com/document/detail/zh/kunpenghpcs/hpckit/instg/KunpengHPCKit_install_012.html)
#### 解压 HPCKit 软件包（HPCKit 版本号根据实际情况调整）
```
tar xvf HPCKit_26.0.RC1_Linux-aarch64.tar.gz
```
#### 安装 HPCKit
```
sh HPCKit_26.0.RC1_Linux-aarch64/install.sh -y --prefix=[HPCKit安装目录]
```

### 3. 安装编译所需依赖

安装 cmake
```
yum install cmake
```

### 4.KuQCD编译安装
KuDNN依赖BiSheng/GCC编译器、HMPI和KUPL，需要预先加载，具体版本号根据HPCKit实际情况调整。这里以26.0.RC1版本HPCKit为例。

#### （1）BiSheng编译KuDNN
##### 加载 module
```
module use [HPCKit安装目录]/HPCKit/latest/modulefiles
```

##### 加载编译器环境变量（编译器版本号根据实际情况调整）
```
module load bisheng/compiler5.1.0.2/bishengmodule
```

##### 加载HMPI环境变量（HMPI版本号根据实际情况调整）
```
module load bisheng/hmpi26.0.RC1/release

export CPATH=[HPCKit安装目录]/HPCKit/latest/hmpi/bisheng/release/hmpi/include:$CPATH

export LIBRARY_PATH=[HPCKit安装目录]/HPCKit/latest/hmpi/hmpi/bisheng/release/hmpi/lib:$LIBRARY_PATH

export LD_LIBRARY_PATH=[HPCKit安装目录]/HPCKit/latest/hmpi/hmpi/bisheng/release/hmpi/lib:$LD_LIBRARY_PATH
```

##### 加载KUPL环境变量（KUPL版本号根据实际情况调整）
```
module load bisheng/kupl26.0.RC1/release
```

##### 编译KuQCD源码
在KuQCD主目录执行相应的编译命令：
```
sh build.sh
```
执行完毕后，会在KuQCD主目录创建build目录并将编译后得到的KuQCD动态库放入该build目录中。

#### （2）GCC编译KuDNN
##### 加载 module
```
module use [HPCKit安装目录]/HPCKit/latest/modulefiles
```

##### 加载编译器环境变量（编译器版本号根据实际情况调整）
```
module load gcc/compiler12.3.1/gccmodule
```

##### 加载HMPI环境变量（HMPI版本号根据实际情况调整）
```
module load gcc/hmpi26.0.RC1/release

export CPATH=[HPCKit安装目录]/HPCKit/latest/hmpi/gcc/release/hmpi/include:$CPATH

export LIBRARY_PATH=[HPCKit安装目录]/HPCKit/latest/hmpi/hmpi/gcc/release/hmpi/lib:$LIBRARY_PATH

export LD_LIBRARY_PATH=[HPCKit安装目录]/HPCKit/latest/hmpi/hmpi/gcc/release/hmpi/lib:$LD_LIBRARY_PATH
```

##### 加载KUPL环境变量（KUPL版本号根据实际情况调整）
```
module load gcc/kupl26.0.RC1/release
```

##### 编译KuQCD源码
在KuQCD主目录执行相应的编译命令：
```
sh build.sh
```
执行完毕后，会在KuQCD主目录创建build目录并将编译后得到的KuQCD动态库放入该build目录中。

### 5. 测试
#### （1）编译运行测试程序(BiSheng)
在已经编译了KuQCD源码后，在KuQCD主目录执行相应的命令：
```
cd test

module load bisheng/hmpi26.0.RC1/release

sh build.sh

sh run.sh
```
#### （2）编译运行测试程序(GCC)
在已经编译了KuQCD源码后，在KuQCD主目录执行相应的命令：
```
cd test

module load gcc/hmpi26.0.RC1/release

sh build.sh

sh run.sh
```

## 📖学习教程

若您已学习编译安装，对本项目有一定认知，并希望深入了解和体验项目，请访问下述详细教程。

[KuQCD开发指南](https://www.hikunpeng.com/document/detail/zh/kunpenghpcs/hpckit/devg/KunpengHPCKit_developer_211.html): 介绍了 KuQCD 中基本数据结构（包括时空格点和格点场），以及求解器相关函数。

## 🔍目录结构

```
├── include                            # 项目工程头文件
├── src                                # 项目源码目录
├── test                               # 项目测试工程目录
├── build.sh                           # 项目工程编译脚本
├── CMakeLists.txt                     # 项目编译配置文件
├── LICENSE                            # LICENSE 文件
└── README.md                          # readme
```

## 🤝联系我们

本项目功能和文档正在持续更新和完善中，建议您关注最新版本。

- **问题反馈**：通过 [【Issues】](https://atomgit.com/kunpengcompute/kuqcd/issues)提交问题。
- **社区互动**：通过 [【鲲鹏论坛（HPC专区）】](https://www.hikunpeng.com/forum/forum-0187135482144798003-1.html)参与交流。
- **技术专栏**：通过 [【鲲鹏社区】](https://www.hikunpeng.com/developer/techArticles) 获取技术文章，如系列化教程、优秀实践等。
