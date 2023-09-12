% Synthetic data set from 
% Minimum-Volume Rank-Deficient Nonnegative Matrix Factorizations, 
% Valentin Leplat, Andersen M.S. Ang, Nicolas Gillis, 2018. 

% Rice REU STAT-DATASCI - Duc Toan Nguyen
% =============
% This code experiments the behaviors of errors of matrix X and matrix W
% with different lambda~ and different noise levels of Min-Vol Algo
% =============

clear all; clc; close all; 

% True generating W, with rank_+(W) = 4 > rank(W) = 3
disp('True W:') 
Wt = [1 0 0 1; 1 0 1 0; 0 1 1 0; 0 1 0 1]'

% Factorization rank 
r = 4; 

% Generate H with n columns
n = 500; 
purity = 0.8; % ~min distance between the columns of W and the data points
alpha = 0.05*ones(r,1);  % Parameter of the Dirichlet distribution
Ht = [sample_dirichlet(alpha,n)']; 
for j = 1 : n
    while max( Ht(:,j) ) > purity
        Ht(:,j) = sample_dirichlet(alpha,1)';
    end
end

%Create matrix Xt with noise = 0
Xt = Wt *Ht;

P_signal = mean(Xt,'all');

%Number of iterations
options.maxiter = 1000;

% For loop of 11 cases with 11 different noises (epsilons)
for j = 1:11

% For each iteration, choose a noise level (epsilon)
if (j==11)
    epsilon=0; %The last case is noise free
else
   epsilon = 10.^(-j); 
end

% Record the noise level to an array eps
eps(j) = epsilon;

% Generate experiment matrix X with noise level epsilon
X = max( 0 , Xt + epsilon*randn(size(Xt)) ); 

SNR(j) = P_signal/epsilon;
SNR_dB(j) = 10*log10(SNR(j));

% Create a list of lambda~ (lambda~ before SNPA initialization)
listK = [1.5 1 0.5 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005 0.00001 0.000005 0.000001 1e-10 1e-11];


W0 = Wt;
minerr=999999;

% Take the loop of all lambda~
for k = 1:length(listK)
    
    % Add lambda~ to the list options
    options.lambda=listK(k);
    
    % Generate the algorithm
    [W,H,e,er1,er2,lambda] = minvolNMF(X,r,options);
    
    trueLambda(k) = options.lambda;

    %Error(Wt,West)
    errorW(k) = compareWs(Wt, W);
    
    %Error(X,Xest)
    errorX(k) = compareWs(Xt,W*H);
end

%==== Result after running with all lambda~ ====
[We_noi(j),indexW] = min(errorW);  %[smallest errorW, index of the smallest errorW] 
[Xe_noi(j),indexX] = min(errorX);  

%We2_noi = errorW with smallest lambda~ 
%(best expected lambda~ for errorX)
We2_noi(j) = errorW(end);         

% bestLambdaW is the best lambda~ for errorW
bestLambdaW(j) = trueLambda(indexW);

%We expect the best lambda~ for errorX is the smallest one.
%If it is not the smallest, we can check the 
%behavior of errorX for each lambda~
if errorX(indexX)<errorX(end)
    bestLambdaX(j) = trueLambda(indexX);

    %Display the difference between the error(X,Xest)
    %with the best lambda and the smallest lambda
    disp("Alert! errorX(end)-errorX(indexX) is "+(errorX(end)-errorX(indexX)));
    str = 'epsilon = '+string(epsilon);
    %Plot the errorX correspond with all lambda~
    figure; 
    plot(log10(listK),log10(errorX),"o-");
    title(str);
    xlabel('log10(lambda~)');
    ylabel('log10(ErrorX)');
    legend('||X-WH||_F^2');
else
    bestLambdaX(j) = trueLambda(end);
end

% ===== Result summary for each iteration =====
disp("Noise level: "+eps(j));
disp("Min Werror with noise "+ eps(j) +" and smallest lambda~ = "+ We2_noi(j));
disp("Min Werror with noise "+ eps(j) +" = "+ We_noi(j));
disp("Best lambda for Werror with noise "+ eps(j) +" = "+bestLambdaW(j));
disp("Min ||X-WH|| with noise "+ eps(j) +" = "+ Xe_noi(j))
disp("Best lambda for Xerror with noise "+ eps(j) +" = "+bestLambdaX(j));
display("==========");
end

%% Plotting

% === Plot the errors with respect to the noise levels (in log scale) ===
figure;
plot(log10(eps),log10(We_noi),'r-+');
hold on;
%plot(log10(eps),log10(We2_noi),'k--o');
plot(log10(eps),log10(Xe_noi),'b--o');
plot(log10(eps),log10(eps),'g-.x');
grid on;
title("$\textbf{Errors}/\textbf{Noise levels}(\sigma) $",'Interpreter','latex');
xlabel("$\log_{10}(\sigma)$",'Interpreter','latex');
ylabel("$\log_{10}(\textbf{errors})$",'Interpreter','latex');
legend('ErrorW','ErrorX','$\sigma$','Interpreter','latex'); 

%% plot best lambda~

% === Plot the best lambda~ for errorX and errorW (in log scale)===
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultLineLineWidth', 2);
figure;
plot(log10(eps),log10(bestLambdaW),'+-');
hold on;
plot(log10(eps),log10(bestLambdaX),'o--');
xlabel("$\log_{10}(\sigma)$",'Interpreter','latex');
ylabel("$\log_{10}(\tilde{\lambda})$",'Interpreter','latex');
grid on;
title('Best $\tilde{\lambda}$','Interpreter','latex');
legend('Best $\tilde{\lambda}_{W}$','Best $\tilde{\lambda}_{X}$','Interpreter','latex');

%%
% === Plot data in 3D space ===
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultLineLineWidth', 2);
figure; 
plot3(X(1,:), X(2,:), X(3,:),'bo'); hold on; 
plot3(Wt(1,:), Wt(2,:), Wt(3,:), 'ro', 'MarkerSize', 12); 
plot3(W(1,:), W(2,:), W(3,:),'mx', 'MarkerSize', 12); 
grid on; 
legend('Data points', 'True W', 'Estimated W'); 

%% 2D plot
% === Plot data in 2D using PCA ===
[U,S,V] = svd(X,'econ');
U1 = U(1:2,:);
X_p = U1*X;
Wt_p = U1*Wt;
W_p = U1*W;

set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultLineLineWidth', 2);
figure; 
scatter(X_p(1,:), X_p(2,:),80,'bo'); hold on; 
scatter(Wt_p(1,:), Wt_p(2,:),100, 'ro',LineWidth=2); 
scatter(W_p(1,:), W_p(2,:),100,'mx',LineWidth=2);  
grid on; 
legend('Data points', 'True W', 'Estimated W'); 



% figure; 
% plot(er1,"o-"); hold on; 
% plot(er2,"+--"); 
% legend('||X-WH||_F^2','logdet(W^TW+ \delta I)'); 
% xlabel('Iterations'); 
% disp('Computed W:')
% 
% %fprintf('Error ||W-Wt||/||Wt|| = %2.2f%%.\n', 100*compareWs( Wt, W ) ); 