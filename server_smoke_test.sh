python python-script/ConfGen.py -i domain.xml || \
	{ echo "usage: testing/domain.xml file in install-direcory legen (über bin und python-script-Ordnern) und dieses Skript von da aus ausführen." ; exit -1; }
mpiexec -n 2 ./bin/SOMA --coord-file domain.h5 --timesteps=10 || exit
mpiexec -n 4 ./bin/SOMA --coord-file domain.h5 --timesteps=10 || exit
mpiexec -n 4 ./bin/SOMA --server-ranks=0,3 --coord-file domain.h5 --timesteps=10 || exit
mpiexec -n 4 ./bin/SOMA --server-ranks=1,2 --N-domains=2 --coord-file domain.h5 --timesteps=10 || exit
mpiexec -n 3 ./bin/SOMA --server-ranks=1 --N-domains=2 --coord-file domain.h5 --timesteps=10 || exit
printf "\n===all tests successful===\n"
