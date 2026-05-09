 # Copyright (c) 2026 Huawei Technologies Co., Ltd. All Rights Reserved.
 #
 # KuQCD is licensed under Mulan PSL v2.
 # You can use this software according to the terms and conditions of the Mulan PSL v2.
 # You may obtain a copy of Mulan PSL v2 at:
 #       http://license.coscl.org.cn/MulanPSL2
 # THIS SOFTWARE IS PROVIDED ON AN "AS IS" BASIS, WITHOUT WARRANTIES OF ANY KIND, EITHER EXPRESS OR IMPLIED,
 # INCLUDING BUT NOT LIMITED TO NON-INFRINGEMENT, MERCHANTABILITY OR FIT FOR A PARTICULAR PURPOSE.
 # See the Mulan PSL v2 for more details.

set(KUQCD_OBJS
        command_line.cpp
        logger.cpp
        kuqcd_utils.cpp
        space_time.cpp
        kuqcd_complex.cpp
        index_utils.cpp
        spinor.cpp
        su3.cpp
        kuqcd_math.cpp
        kuqcd_memory.cpp
        latt_base.cpp
        lattice_field.cpp
        spinor_field.cpp
        gauge_field.cpp
        communication.cpp
        dirac_hisq.cpp
        inverter_cgm.cpp
        kuqcd.cpp
)
# set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Werror -Wall -Wno-instantiation-after-specialization")
# set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror -Wall -Wno-instantiation-after-specialization")
add_library(kuqcd SHARED ${KUQCD_OBJS})
if(USE_SVE)
    target_compile_definitions(kuqcd PRIVATE -DSVEBITS=${SVEBITS})
endif()
if(USE_KUPL)
    target_link_libraries(kuqcd kupl)
endif()

target_link_libraries(kuqcd mpi)
set(HPCKit_VERSION $ENV{HPCKit_VERSION})
set_target_properties(kuqcd PROPERTIES VERSION 1.${HPCKit_VERSION} SOVERSION 1)
