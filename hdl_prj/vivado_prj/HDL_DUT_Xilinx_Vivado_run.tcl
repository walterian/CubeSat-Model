# HDL Coder Downstream Integration Tcl Script

set myTool "Xilinx Vivado 2019.2"
set myProject "HDL_DUT_vivado"
set myProjectFile "HDL_DUT_vivado.xpr"
set myTopLevelEntity "HDL_DUT"

set Family "Virtex7"
set Device "xc7vx485t"
set Package "ffg1761"
set Speed "-2"
set MapTimePathNumber "3"
set MapTimeAdvAnalysis "True"
set PARTimePathNumber "3"
set PARTimeAdvAnalysis "True"

# Open existing project
puts "### Open existing $myTool project hdl_prj\\vivado_prj\\HDL_DUT_vivado.xpr"
open_project ${myProject}

# Running Synthesis
puts "### Running Synthesis in $myTool ..."
if { [llength [get_runs synth_1] ] > 0 } { reset_run synth_1}
launch_runs synth_1 -force
wait_on_run synth_1
puts "### Synthesis Complete."

# Running PostMapTiming
puts "### Running PostMapTiming in $myTool ..."
open_run synth_1
report_timing_summary -delay_type min_max -report_unconstrained -check_timing_verbose -max_paths $MapTimePathNumber -nworst $MapTimePathNumber -name timing_post_map -file timing_post_map.rpt
source -quiet "E:\\MATLAB\\R2019b\\toolbox\\hdlcoder\\hdlcommon\\+downstreamtools\\extractVivadoTiming.tcl"
extracTimingReport ${myTopLevelEntity}_preroute.tvr
puts "### PostMapTiming Complete."

# Close project
puts "### Close $myTool project."
close_project

