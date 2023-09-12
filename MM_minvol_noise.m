% This is the code for 2D ploting best W based on errorW and errorX
%
% Pseudo code:
% +) INPUT: X,r,number of iterations, lambda, delta, 
%        varepsilon (for approximation function)
% +) OUTPUT: W,H
%
% (W0,H0) <- SNPA(X,r)
%
% for i = 1 to number of iteration
%       [W,H] = minvol2(X,r,[W0,H0,lambda,delta])
%       lambda <- 2* sqrt(sum(sum(X-WH))+varepsilon)*lambda  %update lambda
%       W0 <- W
%       H0 <- H
% end for
%
%

%% clear data
clear all; clc; 
%close all; 

rng('default')
s = rng

%% synthetic data

% True generating W, with rank_+(W) = 4 > rank(W) = 3
%disp('True W:') 
Wt = [1 0 0 1; 1 0 1 0; 0 1 1 0; 0 1 0 1]' 
m = 10;
r = 5; % nonnegative rank 
%Wt = rand(m, r);

[m,r] = size(Wt);

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

%% Noiseless X
Xt = Wt*Ht;
normX = norm(Xt,'fro');

P_signal = mean(Xt,'all');


for i = 1:10
noiselevel(i) = 10.^(-i);
 
X = max(0,Xt+noiselevel(i)*randn(size(Xt)));


SNR(i) = P_signal/noiselevel(i);
SNR_dB(i) = 10*log10(SNR(i));


% Run min-vol NMF (default parameters: delta=0.1, maxiter=100)
options = [];
options.lambda = 1;
options.maxiter = 1e3;
epsilon = 0.1;
errorX0=[];
W_F = [];
H_F = [];
normW = norm(Wt,'fro');


for iter = 1:100
    lambda1(iter) = options.lambda;
    [W,H] = minvolNMF2(X,r,options);
    options.lambda = 2*sqrt(sum(sum((X - W*H).^2)) + epsilon)*options.lambda;
    options.W = W;
    options.H = H;
    
    errorX0(iter) = norm(Xt-W*H,'fro')/normX;
    errorW0(iter) = compareWs(W,Wt);

    disp("========Iteration "+iter+"=========");
    disp("ErrorW = "+errorW0(iter));
    disp("ErrorX = "+errorX0(iter));
    disp("Lambda = "+lambda1(iter));

    if iter == 1 
        W_F = W;
        H_F = H;
    else
        if errorX0(iter) >= errorX0(iter-1)
            break;
        else
            W_F = W;
            H_F = H;
        end
    end

end
W0 = W_F;
H0 = H_F;

allerrorX(i) = min(errorX0);
allerrorW(i) = min(errorW0);

disp("=========== SNR ==========");
disp("Noise level: "+noiselevel(i));
disp("SNR = P_signal/noiselevel = "+SNR(i));
disp("SNR_dB = 10log10(SNR) = "+SNR_dB(i));

end



%% SNR
disp("=========== SNR ==========");
disp("Noise level: "+noiselevel(i));
disp("SNR = P_signal/noiselevel = "+SNR(i));
disp("SNR_dB = 10log10(SNR) = "+SNR_dB(i));

%% Plot Errors
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultLineLineWidth', 2);
figure;
plot(log10(noiselevel),log10(allerrorW),'r-+');
hold on;
plot(log10(noiselevel),log10(allerrorX),'b--o');
plot(log10(noiselevel),log10(noiselevel),'g-.*')
grid on;
xlabel("$\log_{10}(\sigma)$",'Interpreter','latex');
ylabel("$\log_{10}(\textbf{errors})$",'Interpreter','latex');
legend("ErrorW", "ErrorX","Noise level");
title("$\textbf{Errors}/\textbf{Noise levels}(\sigma) $",'Interpreter','latex');

%% Plot Lambda

set(0, 'DefaultLineLineWidth', 2);
figure;
plot(log10(lambda1),'*-');
xlabel("iteration");
ylabel("log10(Lambdas)");
title("Lambda/iteration");

%% ===2D plot using PCA===
[U,S,V] = svd(Xt,'econ');
size(U)
size(V)

U1 = U(1:2,:);

X_p = U1*X;
Wt_p = U1*Wt;
W0_p = U1*W0;

set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultLineLineWidth', 2);
figure; 
scatter(X_p(1,:), X_p(2,:),50,'bo'); hold on; 
scatter(Wt_p(1,:), Wt_p(2,:),80, 'ro','LineWidth', 2); 
scatter(W0_p(1,:), W0_p(2,:),80,'mx','LineWidth', 2); 
grid on; 
xlabel("iteration");
ylabel("log10(errors)");
legend('Data points','True W','MinVol W'); 
title("2D plot");

