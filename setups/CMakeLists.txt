
set(setup_scripts
	jureca.sh
	jureca_python.sh
	jureca_cpu.sh
	GAUG_intern.sh
	titan_python.sh
	titan.sh
	hlrn.sh
	strand.sh
	juwels.sh
	juwels_cpu.sh
	juwels_booster.sh
	midway.sh
	midway_cpu.sh)

foreach(SCRIPT ${setup_scripts})
	install(FILES ${SCRIPT} DESTINATION setups)
endforeach(SCRIPT ${setup_scripts})
