SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
LIB_PATH="$SCRIPT_DIR/../build/lib"
LIB_ABS_PATH=$(readlink -f "$LIB_PATH")
INC_PATH="$SCRIPT_DIR/../include"
INC_ABS_PATH=$(readlink -f "$INC_PATH")

export LD_LIBRARY_PATH="$ABS_PATH:${LD_LIBRARY_PATH:-}"
cputype=$(cat /proc/cpuinfo | grep 'CPU part' | awk 'NR==1{print $4}')
    if [[ "$cputype" == "0xd22" ]]; then
        export LIBRARY_PATH="$LIB_ABS_PATH:${LIBRARY_PATH:-}"
        export CPLUS_INCLUDE_PATH="$INC_ABS_PATH:${CPLUS_INCLUDE_PATH:-}"
        mpicxx -lkuqcd test_solver.cpp -o ./test_solver
    fi
