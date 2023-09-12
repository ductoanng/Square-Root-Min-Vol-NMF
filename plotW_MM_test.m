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

%% clear data
clear all; clc; 
%close all; 

rng('default')
s = rng

%% synthetic data

% True generating W, with rank_+(W) = 4 > rank(W) = 3
%disp('True W:') 
%Wt = [1 0 0 1; 1 0 1 0; 0 1 1 0; 0 1 0 1]' 
m = 25;
r = 20; % nonnegative rank 
Wt = rand(m, r);

[m,r] = size(Wt);

% Generate H with n columns
n = 10000; 
purity = 0.8; % ~min distance between the columns of W and the data points
alpha = 0.05*ones(r,1);  % Parameter of the Dirichlet distribution
Ht = [sample_dirichlet(alpha,n)']; 
for j = 1 : n
    while max( Ht(:,j) ) > purity
        Ht(:,j) = sample_dirichlet(alpha,1)';
    end
end

%% Noiseless X
noiselevel = 0;
Xt = Wt*Ht; 
X = max(0,Xt+noiselevel*randn(size(Xt)));
normX = norm(Xt,'fro');

P_signal = mean(Xt,'all');

SNR = P_signal/noiselevel;
SNR_dB = 10*log10(SNR);


% %% Real Data
% load("Sandiego.mat")
% 
% Xt = [];
% %Choosing matrix 4,32,116,128,150 and vectorize it into row vectors and 
% %combine them into a matrix data X with size 5x160000 
% %(each image is 400x400 pixel)
% list = [4,32,116,128,150];
% for i = 1:length(list)
%     M = Sandiego(:,:,list(i));
%     v = M(:);
%     Xt = [Xt;v'];
% end
% 
% X =Xt;
% [m,n] = size(Xt);
% r = 5;
% normX = norm(Xt,'fro');

%% Code



% Run min-vol NMF (default parameters: delta=0.1, maxiter=100)
options = [];
options.lambda = 0.8;
options.maxiter = 1e3;
epsilon = 0.1;
errorX0=[];
W_F = [];
H_F = [];
normW = norm(Wt,'fro');


for iter = 1:100
    options.iteration = 2;
    lambda1(iter) = options.lambda;
    [W,H] = minvolNMF2(X,r,options);
    options.lambda = 2*sqrt(sum(sum((X - W*H).^2)) + epsilon)*options.lambda;
    options.W = W;
    options.H = H;
    
    errorX0(iter) = norm(Xt-W*H,'fro')/normX;
    errorW0(iter) = compareWs(Wt,W);

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

%% SNR
disp("=========== SNR ==========");
disp("Noise level: "+noiselevel);
disp("SNR = P_signal/noiselevel = "+SNR);
disp("SNR_dB = 10log10(SNR) = "+SNR_dB);

%% Plot Errors
set(0, 'DefaultLineLineWidth', 2);
figure;
plot(log10(errorW0),'*-');
hold on;
plot(log10(errorX0),'o-');
grid on;
title("$\textbf{Errors}$/iteration",'Interpreter','latex');
xlabel("iteration",'Interpreter','latex');
ylabel("$\log_{10}(\textbf{errors})$",'Interpreter','latex');
legend('ErrorW','ErrorX','Interpreter','latex'); 

%% Plot Lambda

set(0, 'DefaultLineLineWidth', 2);
figure;
plot(log10(lambda1),'*-');
grid on;
xlabel("iteration",'Interpreter','latex');
ylabel("$log_{10}(\lambda)$",'Interpreter','latex');
title("$\lambda$/iteration",'Interpreter','latex');

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
scatter(Wt_p(1,:), Wt_p(2,:),100, 'ro','LineWidth', 2); 
scatter(W0_p(1,:), W0_p(2,:),100,'mx','LineWidth', 2); 
grid on; 
legend('Data points $X$','$W^{\ast}$ columns','$\hat{W}$ columns','Interpreter','latex'); 
title("2D plot ($\lambda = 1$)",'Interpreter','latex');

% %% Construct H0
% figure;
%  h    = [];
%  h(1) = subplot(4,2,1);
%  h(2) = subplot(4,2,2);
%  h(3) = subplot(4,2,3);
%  h(4) = subplot(4,2,4);
%  h(5) = subplot(4,2,5);
%  % h(6) = subplot(4,2,6);
%  % h(7) = subplot(4,2,7);
%  % h(8) = subplot(4,2,8);
% 
% for j=1:r
%     imagesc(reshape(H0(j,:),400,400),'Parent',h(j));
% end