clean:
	rm -f *.nc *.nc.* MOM_parameter_doc.* mom6.standalone.o* ocean_solo.res exitcode available_diags.000000 ocean.stats time_stamp.out logfile.000000.out CPU_stats

rose-clean: clean
	# remove rose generated files
	rm -rf prose_logs prose_workspace *postprocessed* *.rmod rose_*

load: rose-clean
	cp -r ../../MOM6_full_fp_dataflow_graph_4/* .

reset: rose-clean
	rm ../../build/* -rf