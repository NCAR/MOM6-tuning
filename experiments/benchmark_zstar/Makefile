TARGET_DIRECTORY := ../../src/MOM6/src/core

clean:
	rm -f *.nc *.nc.* MOM_parameter_doc.* mom6.standalone.o* ocean_solo.res exitcode available_diags.000000 ocean.stats time_stamp.out logfile.000000.out CPU_stats

rose-clean: clean
	# remove rose generated files
	rm -rf prec_logs prec_workspace *postprocessed* *.rmod 

reset: rose-clean
	rm ../../build/* -rf