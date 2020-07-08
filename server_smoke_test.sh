python python-script/ConfGen.py -i server.xml || \
	{ echo "usage: testing/server.xml file in install-direcory legen (über bin und python-script-Ordnern) und dieses Skript von da aus ausführen." ; exit -1; }
mpiexec --tag-output -n 2 ./bin/SOMA --coord-file server.h5 --ana-file server_ana.h5 --timesteps=10 || exit
mpiexec -n 4 ./bin/SOMA --coord-file server.h5 --ana-file server_ana.h5 --timesteps=500 || exit
mpiexec -n 4 ./bin/SOMA --server-ranks=0,3 --ana-file server_ana.h5 --coord-file server.h5 --timesteps=500 || exit
mpiexec -n 4 ./bin/SOMA --server-ranks=1,2 --ana-file server_ana.h5 --N-domains=2 --coord-file server.h5 --timesteps=500 || exit
mpiexec -n 3 ./bin/SOMA --server-ranks=1 --N-domains=2 --ana-file server_ana.h5 --coord-file server.h5 --timesteps=500 || exit
printf "\n===all tests successful===\n"
