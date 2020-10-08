TARGET_DIRECTORY := ../../src/MOM6/src/core
TARGET_FILE := MOM_CoriolisAdv.F90

clean:
	rm -f *.nc *.nc.* MOM_parameter_doc.* mom6.standalone.o* ocean_solo.res exitcode available_diags.000000 ocean.stats time_stamp.out logfile.000000.out CPU_stats

rose-clean: clean
	# remove rose generated files
	rm -rf $(TARGET_FILE) prec_logs prec_workspace *postprocessed* *.rmod initial_config.txt

	# remove the preprocessed target file in the target directory
	rm -f $(TARGET_DIRECTORY)/preprocessed_$(TARGET_FILE)

reset: rose-clean
	rm ../../build/* -rf