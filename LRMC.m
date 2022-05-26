function [Q]=LRMC(P, Omega)
% this function completes a low rank matrix P with known entries set Omega
[n1,n2] = size(P);
maxiters=30; 
% maximum number of iterations
epsilon = 1e-5;                 
tau = 1.5;
%% Calculating noise variance of patch P
flag_second_subset = 1;
% set flag_secon_subset = 1 if want to consider second set of noisy pixels
% in misssing set else 0
std_dev_bar = zeros(n1, 1);
for i=1:n1
    P_ith_row = P(i,:);
    std_dev_bar(i) = std(P_ith_row(Omega(i,:)));
    if flag_second_subset
        mu_bar = mean(P_ith_row(Omega(i,:)));            
        Omega(i,:) = Omega(i,:) & (P_ith_row <= mu_bar + 2 * std_dev_bar(i)) & (P_ith_row >= mu_bar - 2 * std_dev_bar(i));
    end
end
std_dev_hat = sum(std_dev_bar,"all")/numel(std_dev_bar);
%% LRMC algoritm
p = sum(Omega, 'all')/numel(Omega);
mu = (sqrt(n1) + sqrt(n2)) * sqrt(p) * std_dev_hat;
Q = zeros(size(P));
iter = 1;
Q_next = Q;
while(sqrt(sum((Q_next-Q).^2,"all"))<=epsilon&&iter<=maxiters)
%     while iterations <= maximum iterations or mse is less than epsilon
    P_proj = (Q-P).*Omega;
    R = Q-tau*P_proj;  
    [U, S, V] = svd(R, 'econ');
    Q_next = U * max(S - tau*mu, 0) * V.'; 
%     Soft shrinkage
    Q = Q_next; 
    iter = iter + 1;
end   
end
