python python-script/ConfGen.py -i coord.xml || \
	{ echo "usage: coord.xml file in install-direcory legen (über bin und python-script-Ordnern) und dieses Skript von da aus ausführen." ; exit -1; }
mpiexec -n 2 ./bin/SOMA --coord-file coord.h5 --timesteps=10
mpiexec -n 4 ./bin/SOMA --coord-file coord.h5 --timesteps=10
mpiexec -n 4 ./bin/SOMA --N-servers=2 --coord-file coord.h5 --timesteps=10
mpiexec -n 4 ./bin/SOMA --N-servers=2 --N-domains=2 --domain-buffer=4 --coord-file coord.h5 --timesteps=10
mpiexec -n 3 ./bin/SOMA --N-domains=2 --domain-buffer=4 --coord-file coord.h5 --timesteps=10
