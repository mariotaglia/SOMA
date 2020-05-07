clear
make || exit 1
mpirun  -np 2 ./runtests --gtest_filter=*t2* || exit 1
mpirun --oversubscribe -np 4 ./runtests || exit 1
mpirun --oversubscribe -np 7 ./runtests || exit 1
mpirun --oversubscribe -np 16 ./runtests || exit 1
echo "All tests succeeded"
