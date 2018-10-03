#/bin/sh
module add R
echo "simulations 1,2"
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data") &
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_2") &
wait
echo "simulations 3,4"
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_3") &
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_4") &
wait
echo "simulation 5"
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_5") &
echo "simulated_data_6"
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_6") &
wait
echo "simulated_data_7,8"
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_7") &
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_8") &
wait
echo "simulated_data_9,10"
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_9") &
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_10") &
wait
echo "simulated_data_11,12"
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_11") &
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_12") &
wait
echo "simulated_data_13,14"
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_13") &
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_14") &
wait
echo "simulated_data_15"
(Rscript ./simulations_run_post_analysisU.r "./simulationsU/simulated_data_15") &
wait
echo "running plot"
Rscript simulations_run_plot.r "./simulationsU"
