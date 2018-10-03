
function simulations_run_analysis(simdata_dir)

% The analysis to run (mark 0 to ignore):
run_estimate_params = 0;
run_Refactor_fix_p = 0;
run_Refactor_fix_tau = 0;
run_robustness_t = 0;
run_robustness_K = 0;
run_variance_R =1
addpath './'

% Higher values of tau makes it easier for a standard PCA to find R.
tau_values = 0.01:0.02:0.1;
% Higher values of p makes it easier for a standard PCA to find R.
p_values = 0.1:0.2:1;
p_fixed = 0.15;
% Number of individuals to simulate in each dataset.
n = 2000;
% Noise of the technology.
sigma = 0.01;

sorted_cells = './m_sorted_103638.mat';
estimates_dir = './parameter_estimates/';
%simdata_dir = './simulations/simulated_data/';
results_dir = './simulations/';

% Parameters for Refactor.
MAX_COMPONENTS = 20;
T = 500;
K = 5;

% First, estimate tau from the sorted cells data, under the assumption that
% all sites are DMR sites.
load(sorted_cells);
M = m';
clear m;
means = zeros(5,size(M,2));
for k = 1:5
    means(k,:) = mean(M((k-1)*5+1:k*5,:));
end
tau_est = (sum(sum((means-repmat(mean(means),5,1)).^2)) / 518190).^0.5;
tau_fixed = tau_est - mod(tau_est,0.01);


%% (1) Estimate parameters from Reinius reference.

if (run_estimate_params)
    
    display('run_estimate_params...')
    for tau = tau_values
        display(tau)
        outfile = [estimates_dir 'parameters_estimated_from_sorted_tau_' num2str(tau) '.mat'];
        estimate_parameters(tau, outfile);
    end
    
end


%% (2) Fix p and simulate data using different tau values.

% According to Jaffe 2014, 15% of the sites in the illumina 450K chip
% passed their test for begin DMR sites (after Bonferroni correction).

if (run_Refactor_fix_p)
    
    display('run_Refactor_fix_p...')
    
    resultsfile = [results_dir 'sim_results_fixed_p_' num2str(p_fixed) '.mat'];
    
    refactor_results_fixed_p = cell(1,length(tau_values));
    pca_results_fixed_p = cell(1,length(tau_values));
    
    for tau = tau_values
        display(tau)
        estimates_file = [estimates_dir 'parameters_estimated_from_sorted_tau_' num2str(tau) '.mat'];
        outfile = [simdata_dir 'sim_data_n_' num2str(n) '_p_' num2str(p_fixed) '_tau_' num2str(tau) '.mat'];
        simulate_data_hierarchical_model(n,p_fixed,sigma,estimates_file,outfile,1000);
    end
end



%% (3) Fix tau and simulate data using different values of p.


if (run_Refactor_fix_tau)
    
    display(['run_Refactor_fix_tau with tau: ' num2str(tau_fixed)]);
    resultsfile = [results_dir 'sim_results_fixed_tau_' num2str(tau_fixed) '.mat'];
    refactor_results_fixed_tau = cell(1,length(p_values));
    pca_results_fixed_tau = cell(1,length(p_values));
    
    estimates_file = [estimates_dir 'parameters_estimated_from_sorted_tau_' num2str(tau_fixed) '.mat'];
    for p = p_values
        display(p)
        outfile = [simdata_dir 'sim_data_n_' num2str(n) '_p_' num2str(p) '_tau_' num2str(tau_fixed) '.mat'];
        simulate_data_hierarchical_model(n,p,sigma,estimates_file,outfile,1000);
    end

end

%added by daniel
if(run_variance_R)
estimates_file = [estimates_dir 'parameters_estimated_from_sorted_tau_' num2str(tau_fixed) '.mat'];
    for sp = [0.001,0.1,1,10,1000]
        display(sp)
        outfile = [simdata_dir 'sim_data_n_' num2str(n) '_p_' num2str(p_fixed) '_tau_' num2str(tau_fixed) '_sp_' num2str(sp) '.mat'];
        simulate_data_hierarchical_model(n,p_fixed,sigma,estimates_file,outfile,sp);
    end
end


%% (4) Robustness of Refactor to the selection of t

if (run_robustness_t)
    
    %estimates_file = [estimates_dir 'parameters_estimated_from_sorted_tau_' num2str(tau_fixed) '.mat'];
    %outfile = [simdata_dir 'sim_data_n_' num2str(n) '_p_' num2str(p_fixed) '_tau_' num2str(tau_fixed) '.mat'];
    %simulate_data_hierarchical_model(n,p,sigma,estimates_file,outfile);
    
    display('run_robustness_t...')
    
    resultsfile = [results_dir 'sim_results_robustness_t_tau_' num2str(tau_fixed) '_p_' num2str(p_fixed) '.mat'];
    simfile = [simdata_dir 'sim_data_n_' num2str(n) '_p_' num2str(p_fixed) '_tau_' num2str(tau_fixed) '.mat'];
    load(simfile);
    T_values = [100:100:1000 2000:1000:10000 15000:5000:30000];
    [R_est,~] = Refactor(O,K,T_values,MAX_COMPONENTS);
    [~,score] = pca(zscore(O'));
    
    pca_results_robust_t = zeros(MAX_COMPONENTS,5);
    refactor_results_robust_t = cell(length(T_values),1);
    for c = 1:MAX_COMPONENTS
        for k = 1:5
            pca_results_robust_t(c,k) = lin_fit(R(:,k),score(:,1:c));
        end
    end
        
    for t = 1:length(T_values)
        corrs_refactor = zeros(MAX_COMPONENTS,5);
        for c = 1:MAX_COMPONENTS
            for k = 1:5
                corrs_refactor(c,k) = lin_fit(R(:,k),R_est{t}(:,1:c));
            end
        end
        refactor_results_robust_t{t} = corrs_refactor;
    end
    save(resultsfile,'refactor_results_robust_t','pca_results_robust_t','T_values');
    
end



%% (5) Robustness of Refactor to the selection of K

if (run_robustness_K)
    
    %estimates_file = [estimates_dir 'parameters_estimated_from_sorted_tau_' num2str(tau_fixed) '.mat'];
    %outfile = [simdata_dir 'sim_data_n_' num2str(n) '_p_' num2str(p_fixed) '_tau_' num2str(tau_fixed) '.mat'];
    %simulate_data_hierarchical_model(n,p,sigma,estimates_file,outfile);
    
    display('run_robustness_K...')
    
    resultsfile = [results_dir 'sim_results_robustness_K_tau_' num2str(tau_fixed) '_p_' num2str(p_fixed) '.mat'];
    simfile = [simdata_dir 'sim_data_n_' num2str(n) '_p_' num2str(p_fixed) '_tau_' num2str(tau_fixed) '.mat'];
    load(simfile);
    K_values = 3:10;
    [~,score] = pca(zscore(O'));
    pca_results_robust_K = zeros(MAX_COMPONENTS,5);
    refactor_results_robust_K = cell(length(K_values),1);
    for c = 1:MAX_COMPONENTS
        for k = 1:5
            pca_results_robust_K(c,k) = lin_fit(R(:,k),score(:,1:c));
        end
    end
    
    counter = 1;
    for k_value = K_values
        [R_est,~] = Refactor(O,k_value,T,MAX_COMPONENTS);
        
        corrs_refactor = zeros(MAX_COMPONENTS,5);
        for c = 1:MAX_COMPONENTS
            for k = 1:5
                corrs_refactor(c,k) = lin_fit(R(:,k),R_est{1}(:,1:c));
            end
        end
        refactor_results_robust_K{counter} = corrs_refactor;
        counter = counter + 1;
    end
    
    save(resultsfile,'refactor_results_robust_K','pca_results_robust_K','K_values');
    
end



end
