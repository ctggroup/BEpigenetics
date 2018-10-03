%% simulate_data simulated DNA methylation data coming from a mixture of cell types.

% Input:
% (1) n is the number of individuals to simulate. 
% (2) p is the proportion of DMR sites (out of the total number of sites).
% (3) sigma is standard deviation of the noise arising from the technology
% (assumed to be normally distributed with mean 0).
% (4) estimates_file is a path to a file containing estimates for the parameter of the model (given a specific tau value).

% Output (saved into a file):
% (1) O is a n individuals by m sites DNA methylation matrix.
% (2) R is the real cell type proportions used for generating O.
% (3) dmrs is a vector with the indices of the dr sites in O.

% For DMR sites the function simulates data following the model:
% O_{ij}~N(sum_k(M_{jk}R_{ki}),(sigma_i^2)*(sum_k{R_{ik}})).
% If i is not a DMR site then M_{i1}=...=M_{iK}.
% The DMR sites are randomly selected.

% Notes:
% Make sure fastfit and lightspped packages arein the path - for estimating
% and sampling from Dirichete distribution.
% Replace the original implementation of randgamma (randgamma.m file) with
% funciton: x = randg(a);

function simulate_data_hierarchical_model(n,p,sigma,estimates_file,outfile,SAMPLE_CONFIDENCE)

addpath 'fastfit'
addpath 'lightspeed'


load(estimates_file);
m = length(lrt);
K = 5;

%% (1) Select DMR sites


% Randomly select sites as DMR sites
sorted_sites = randperm(m);
dmrs = sorted_sites(length(sorted_sites)-round(p*m)+1:length(sorted_sites));
non_dmrs = setdiff(1:m,dmrs);
mus = zeros(m,5);
sigs = zeros(m,1);
mus(dmrs,:) = mus_H1(dmrs,:);
mus(non_dmrs,:) = repmat(mus_H0(non_dmrs)',1,5);
sigs(dmrs,:) = sigs_H1(dmrs,:);
sigs(non_dmrs,:) = sigs_H1(non_dmrs,:);
sigs = repmat(sigs,1,5);


%% (2) Simulate cell type proportions

% counts for learning the Dirichlete distribution - represends the confidence in the estimated of the cell proportions.
LEARN_CONFIDENCE = 1000;
% counts for sampling form a Dirichlete distribution - represends the
% variance level of the sample.
%SAMPLE_CONFIDENCE = 1000;
% Fit a Dirichlete distribution for R using Houseman estimates of 5 cell types for the RA dataset.
h = dlmread('h4mat.txt');
a = polya_fit(round(h*LEARN_CONFIDENCE));
R = zeros(n,size(h,2));
for j = 1:n
    R(j,:) = polya_sample(a,SAMPLE_CONFIDENCE);
    R(j,:) = R(j,:) ./ sum(R(j,:),2);
end


%% (3) Generate DNA methylation levels

% Generate an ovserved methylation matrix.
% Make sure no values outside [0,1] are generated in any cell type.
O = zeros(n,m);
counter = 1;
%R_squared_norm = sum(R.^2,2);
r = zeros(K,m);
for j = (1:n) 
    r = normrnd(mus',sigs',K,m);
    r(r < 0) = 0;
    r(r > 1) = 1;
    O(j,:) = R(j,:)*r;
    % Display progress.
    counter = counter + 1;
    if (mod(counter,100) == 0)
        display(counter)
    end
end
% Add noise at level sigma to simulate the variance arising form the technology.
M = O;
O = O + normrnd(0,sigma,size(O,1),size(O,2));
indices = find(O < 0);
O(O < 0) = zeros(size(indices));
indices = find(O > 1);
O(O > 1) = ones(size(indices));

O = O';

outfile
save(outfile ,'M','O','R','dmrs','p','n','estimates_file','sigma');

end
