# KuQCD

## Release Notes

- [2026/04] The Kunpeng Unified Lattice Quantum Chromodynamics Library (KuQCD) project went live, introducing a Multishift CG solver optimized for HISQ dslash operations.

## Overview

KuQCD is a C/C++ acceleration library for lattice QCD (LQCD), purpose-built for the Kunpeng platform. It currently supports the Kunpeng 920 new model and depends on KUPL. Supporting both BiSheng Compiler and GCC, KuQCD offers essential LQCD modules, including inverters and dslash/solver algorithms specifically tuned to maximize the hardware potential of Kunpeng processors through high-performance APIs.

## Version Mapping

- Operating platform
    - Kunpeng 920 Pro
- System specifications
    - openEuler 22.03 LTS SP4 AArch64

## Build and Installation

### 1. [Obtain the HPCKit software package](https://www.hikunpeng.com/document/detail/zh/kunpenghpcs/hpckit/instg/KunpengHPCKit_install_007.html).

> https://www.hikunpeng.com/developer/hpc/hpckit-download

### 2. [Install HPCKit](https://www.hikunpeng.com/document/detail/zh/kunpenghpcs/hpckit/instg/KunpengHPCKit_install_012.html).

#### Extract the HPCKit software package (replace the version number in the example with your actual version).

```shell
tar xvf HPCKit_26.2.0_Linux-aarch64.tar.gz
```

#### Install HPCKit.

```shell
sh HPCKit_26.2.0_Linux-aarch64/install.sh -y --prefix=[HPCKit_installation_directory]
```

### 3. [Set environment variables](https://www.hikunpeng.com/document/detail/zh/kunpenghpcs/hpckit/instg/KunpengHPCKit_install_014.html).

#### Load the module.

```shell
module use [HPCKit_installation_directory]/HPCKit/latest/modulefiles
```

#### Load the environment variables of the compiler.

Identify your preferred compiler (GCC or BiSheng Compiler) and execute the loading command via the terminal.

- GCC (Replace the version number in the example with your actual version.)

```shell
module load gcc/compiler12.3.1/gccmodule
```

- BiSheng Compiler (Replace the version number in the example with your actual version.)

```shell
module load bisheng/compiler5.1.0.2/bishengmodule
```

#### Load the HMPI environment variables.

- GCC (Replace the HMPI version number in the example with your actual version.)

```shell
module load gcc/hmpi26.2.0/release

export CPATH=[HPCKit_installation_directory]/HPCKit/latest/hmpi/gcc/release/hmpi/include:$CPATH

export LIBRARY_PATH=[HPCKit_installation_directory]/HPCKit/latest/hmpi/hmpi/gcc/release/hmpi/lib:$LIBRARY_PATH

export LD_LIBRARY_PATH=[HPCKit_installation_directory]/HPCKit/latest/hmpi/hmpi/gcc/release/hmpi/lib:$LD_LIBRARY_PATH
```

- BiSheng Compiler (Replace the HMPI version number in the example with your actual version.)

```shell
module load bisheng/hmpi26.2.0/release

export CPATH=[HPCKit_installation_directory]/HPCKit/latest/hmpi/bisheng/release/hmpi/include:$CPATH

export LIBRARY_PATH=[HPCKit_installation_directory]/HPCKit/latest/hmpi/hmpi/bisheng/release/hmpi/lib:$LIBRARY_PATH

export LD_LIBRARY_PATH=[HPCKit_installation_directory]/HPCKit/latest/hmpi/hmpi/bisheng/release/hmpi/lib:$LD_LIBRARY_PATH
```

#### Load the KUPL environment variables.

- GCC (Replace the KUPL version number in the example with your actual version.)

```shell
module load gcc/kupl26.2.0/release
```

- BiSheng Compiler (Replace the KUPL version number in the example with your actual version.)

```shell
module load bisheng/kupl26.2.0/release
```

### 4. Install the dependencies required for build.

Install CMake.

```shell
yum install cmake
```

### 5. Start the build process.

#### Build KuQCD from source.

Run the following build command on the terminal:

```shell
sh build.sh
```

#### Build and run the test program (GCC).

After building the source, run the following commands in the KuQCD main directory:

```shell
cd test

module unload bisheng/hmpi26.2.0/release

module load gcc/hmpi26.2.0/release

sh build.sh

sh run.sh
```

#### Build and run the test program (BiSheng Compiler).

After building the source, run the following commands in the KuQCD main directory:

```shell
cd test

module unload gcc/hmpi26.2.0/release

module load bisheng/hmpi26.2.0/release

sh build.sh

sh run.sh
```

## Tutorials

If you are familiar with the build and installation process and would like to gain a deeper understanding of the project, please visit the following detailed tutorials:

[KuQCD Developer Guide](https://www.hikunpeng.com/document/detail/zh/kunpenghpcs/hpckit/devg/KunpengHPCKit_developer_211.html): Explore foundational data structures, including spacetime lattices and lattice fields, and master solver functions.

## Directory Structure

```shell
├── include                            # Project header files
├── src                                # Project source directory
├── test                               # Project test directory
├── build.sh                           # Project build script
├── CMakeLists.txt                     # Project build configuration file
├── LICENSE                            # License file
└── README.md                          # Readme
```

## Contact Us

Features and documentation are updated regularly. Please follow the latest version for the most up-to-date information.

- **Issue feedback**: Submit queries or report bugs via [Issues](https://atomgit.com/kunpengcompute/kuqcd/issues).
- **Community interaction**: Join discussions and share ideas via [Discussions](https://atomgit.com/kunpengcompute/kuqcd/discussions).
- **Technical columns**: Access in-depth technical articles, serialized tutorials, and best practices through the [Kunpeng Community](https://www.hikunpeng.com/developer/techArticles).
