#/bin/sh
module add R
echo "simulations 1,2"
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data") &
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_2") &
echo "simulations 3,4"
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_3") &
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_4") &
echo "simulation 5"
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_5") &
echo "simulated_data_6"
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_6") &
echo "simulated_data_7,8"
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_7") &
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_8") &
echo "simulated_data_9,10"
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_9") &
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_10") &
echo "simulated_data_11,12"
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_11") &
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_12") &
echo "simulated_data_13,14"
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_13") &
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_14") &
echo "simulated_data_15"
(Rscript ./simulations_run_bayes_no_sparse.R "../simulations_no_sparse/simulated_data_15") &
