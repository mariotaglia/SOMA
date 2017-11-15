
set(setup_scripts
	jureca.sh
	jureca_python.sh
	GAUG_intern.sh
	titan_python.sh
	titan.sh )

foreach(SCRIPT ${setup_scripts})
	install(FILES ${SCRIPT} DESTINATION setups)
endforeach(SCRIPT ${setup_scripts})
