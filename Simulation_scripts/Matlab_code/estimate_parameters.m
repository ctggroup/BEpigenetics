% estimate_parameters function uses the sorted cells methylation levels
% collected by Reinius et al. for estimating the parameters of the
% following hierarchical model assumed for methylation sites.
% For a site j, the methylation levels in cell type k are coming from
% N(mu_jk,sigma_j^2), where mu_jk comes from N(mu_j,tau^2).
%
% The function get one parameter:
% A value for tau - as tau gets higher, the DMR sites distinct better
% between different cell types, thus making it easier for a standard PCA to
% unveil the cell composition.
%
% The function saves a file with the variables:
% mus_H0,mus_H1: mean values for each methylation site at each cell type under the
% null model (not DMR) and under the DMR model.
% sigs_H0,sigs_H1: standard deviation for each methylation site under the
% null model (not DMR) and under the DMR model.
% lrt: log likelihood ratio statistic for each site.

function estimate_parameters(tau,outfile)

display('Start estimate_parameters...')

% The following containes 6 consecutive samples for each cell type (for
% each site).
%load '../Resources/Sorted/m_sorted_103638.mat';
%load '/home/nasheran/eliorrah/backup/Refactor/simulations/sorted_meth/m_sorted_103638.mat';
load './m_sorted_103638.mat'

M = m;
clear m;

m = size(M,1);
J = 6;
K = 5;
% replace nan values (-1) by the mean of that row.
for i = 1:m
    f = find(M(i,:) < 0);
    if (sum(f) > 0)
        M(i,f) = repmat(mean(M(i,:)),1,length(f));
    end
end
M = M';

% Estimate the parameters under the null model.
mus_H0 = mean(M);
sigs_H0 = std(M);

% Estimate the parameters under the alternative model.
mus_H1 = zeros(m,5);
sigs_H1 = zeros(m,1);
warning ('off','all');
counter = 1;
lrt = zeros(m,1);
n = 6;
%tic()
for j = 1:m
    Y = reshape(M(:,j),6,5);
    % Format a fsolve problem, where the input is a vector x, where x(1:5):
    % mu_jk, x(6): sigma_j
    mu_j = sum(sum(Y))/30;
    s1 = sum(Y(:,1));
    s2 = sum(Y(:,2));
    s3 = sum(Y(:,3));
    s4 = sum(Y(:,4));
    s5 = sum(Y(:,5));
    ss1 = sum(Y(:,1).^2);
    ss2 = sum(Y(:,2).^2);
    ss3 = sum(Y(:,3).^2);
    ss4 = sum(Y(:,4).^2);
    ss5 = sum(Y(:,5).^2);
    myfun = @(x) [x(1)-( ((tau^2)*s1+(x(6)^2)*mu_j)/(n*(tau^2)+x(6)^2) ) ;
        x(2)-( ((tau^2)*s2+(x(6)^2)*mu_j)/(n*(tau^2)+x(6)^2)) ;
        x(3)-( ((tau^2)*s3+(x(6)^2)*mu_j)/(n*(tau^2)+x(6)^2)) ;
        x(4)-( ((tau^2)*s4+(x(6)^2)*mu_j)/(n*(tau^2)+x(6)^2)) ;
        x(5)-( ((tau^2)*s5+(x(6)^2)*mu_j)/(n*(tau^2)+x(6)^2)) ;
        x(6)^2-( ( n*x(1)^2-2*x(1)*s1+ss1 + n*x(2)^2-2*x(2)*s2+ss2 + n*x(3)^2-2*x(3)*s3+ss3 + n*x(4)^2-2*x(4)*s4+ss4 + n*x(5)^2-2*x(5)*s5+ss5 )/(5*n) );
        tau^2 - sum((x(1:5)-mu_j).^2)/5;];
    x0 = [mean(Y(:,1));mean(Y(:,2));mean(Y(:,3));mean(Y(:,4));mean(Y(:,5));std(Y(:))];
    %options = optimoptions('fsolve','Display','iter');
    options = optimset('Display','off');
    [x,fval] = fsolve(myfun,x0,options);
    x = abs(x)';
    mus_H1(j,:) = x(1:5);
    sigs_H1(j) = x(6);
    % Calculate the log likelihood ratio of the two models.
    l0 = -30*log(sigs_H0(j))-sum((M(:,j)-mus_H0(j)).^2)/(2*(sigs_H0(j)^2));
    l1 = -30*log(sigs_H1(j)) - 2.5*log(2*pi*(tau^2)) - sum(sum((Y-repmat(mus_H1(j,:),n,1)).^2))/(2*(sigs_H1(j)^2)) - sum((mus_H1(j,:)-mu_j).^2)/(2*(tau^2));
    lrt(j) = l1-l0;
    
    if (mod(counter,10000) == 0)
        display(counter)
    end
    counter = counter + 1;
end
%t = toc()

save(outfile,'mus_H1','sigs_H1','mus_H0','sigs_H0','tau','lrt');



% j = 1;
% 
% V_est = mean((M(:,j) - mean(M(:,j))).^2);
% sigma = sqrt(V_est-tau^2);
% %mu = ((tau^2)*sum(sum(Y))) / ( 5*(n*(tau^2) + sqrt(V_est-tau^2) - sigma^2) );
% %mu = sum(sum(Y)) / (5*(n-(sigma/tau)^2+sigma^2))
% %mu = sum(sum(Y))/30;
% k = 1;


%% % (2)Using fsolve
% for j = 1:m
%     j
%     Y = reshape(M(:,j),6,5);
%     % Format a fsolve problem, where the input is a vector x, where x(1:5):
%     % mu_jk, x(6): sigma_j, x(7): mu_j
%     n = 6;
%     s1 = sum(Y(:,1));
%     s2 = sum(Y(:,2));
%     s3 = sum(Y(:,3));
%     s4 = sum(Y(:,4));
%     s5 = sum(Y(:,5));
%     ss1 = sum(Y(:,1).^2);
%     ss2 = sum(Y(:,2).^2);
%     ss3 = sum(Y(:,3).^2);
%     ss4 = sum(Y(:,4).^2);
%     ss5 = sum(Y(:,5).^2);
%     myfun = @(x) [x(1)-( ((tau^2)*s1+(x(6)^2)*x(7))/(n*(tau^2)+x(6)^2) ) ;
%         x(2)-( ((tau^2)*s2+(x(6)^2)*x(7))/(n*(tau^2)+x(6)^2) ) ;
%         x(3)-( ((tau^2)*s3+(x(6)^2)*x(7))/(n*(tau^2)+x(6)^2) ) ;
%         x(4)-( ((tau^2)*s4+(x(6)^2)*x(7))/(n*(tau^2)+x(6)^2) ) ;
%         x(5)-( ((tau^2)*s5+(x(6)^2)*x(7))/(n*(tau^2)+x(6)^2) ) ;
%         x(6)^2-( ( n*x(1)^2-2*x(1)*s1+ss1 + n*x(2)^2-2*x(2)*s2+ss2 + n*x(3)^2-2*x(3)*s3+ss3 + n*x(4)^2-2*x(4)*s4+ss4 + n*x(5)^2-2*x(5)*s5+ss5 )/(5*n) );
%         x(7)-((x(1)+x(2)+x(3)+x(4)+x(5))/5);];
%     x0 = [mean(Y(:,1));mean(Y(:,2));mean(Y(:,3));mean(Y(:,4));mean(Y(:,5));mean(mean(Y));0.05];
%     %options = optimoptions('fsolve','Display','iter');
%     options = optimoptions('fsolve');
%     %options = optimoptions('fsolve','MaxFunEvals',50000,'TolFun',10^-10,'MaxIter',50000);
%     [x,fval] = fsolve(myfun,x0,options);
%     x = abs(x)';
%     mus_H1(j,:) = x(1:5);
%     sigs_H1(j) = x(6);
% end


% 
% 
% %% An attempt of an analytical solution using the low of total variance (for calculating an
% % estimate for sigma).
% n = 6;
% lrt = zeros(m,1);
% for j = 1:m
%     Y = reshape(M(:,j),6,5);
%     V_est = mean((M(:,j) - mean(M(:,j))).^2);
%     % From the low of total variation we get:
%     sigma_j = sqrt(V_est-tau^2);
%     for k = 1:5
%         mus_H1(j,k) = ( (tau^2)*(sum(Y(:,k)))-(sigma_j^2)*(sum(sum(Y)))/30 ) / (n*(tau^2) - sigma_j^2);
%     end
%     sigs_H1(j) = sigma_j;
%        
%     % Find the log likelihood ratio.
%     mu_j = mean(mus_H1(j,:));
%     
%     mu_j_H0 = mean(M(:,j));
%     sigma_j_H0 = (sum((M(:,j)-mu_j_H0).^2) / 30)^0.5;
%     
%     l0 = -30*log(sigma_j_H0)-sum((M(:,j)-mu_j_H0).^2)/(2*(sigma_j_H0^2));
%     l1 = -30*log(sigs_H1(j)) - 2.5*log(2*pi*(tau^2)) - sum(sum((Y-repmat(mus_H1(j,:),n,1)).^2))/(2*(sigs_H1(j)^2)) - sum((mus_H1(j,:)-mu_j).^2)/(2*(tau^2));
%     lrt(j) = l1-l0;
%     
% 
% end

end
