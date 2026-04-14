SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)          
LIB_PATH="$SCRIPT_DIR/../build/lib"
LIB_ABS_PATH=$(readlink -f "$LIB_PATH")
export LD_LIBRARY_PATH="$LIB_ABS_PATH:${LD_LIBRARY_PATH:-}"
output=$(./test_solver < input)
echo "start test kuqcd_solve_multishift"
success=$(echo "$output" | grep -oi 'PASSED' | wc -l)
if [ $success -eq 2 ]; then
    echo "end test kuqcd_solve_multishift PASSED"
    echo "1 tests ran, 1 tests PASSED, 0 tests FAILED"
else
    echo "end test kuqcd_solve_multishift FAILED"
    echo "1 tests ran, 1 tests FAILED, 0 tests PASSED"
fi