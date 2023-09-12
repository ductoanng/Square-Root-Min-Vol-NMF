% TEST between TITANized Minvol and MM Minvol
% DATASET: Urban
%==========================================================
%% clear data
clear all;
clc;


%% load data
dataset = "urban";

if dataset == "urban"
    load('Urban.mat')
    X = A';
    clear A
    r = 6;
    %maxtime = 60;
end

%% TITANized Minvol

% Initialization 
options = [];
options.lambda=0.001;
options.display = true;
options.inneriter = 10;
options.delta_iter = 1e-2;
options.save = true;
n_samples = 300; % If options.save == true, number of points in the .dat file
n_runs = 20; % Number of runs tested per dataset

k = 0;
nX = norm(X,'fro');
[m,n] = size(X);

W = rand(m,r);
H = FGMfcnls(X,W,[],100);   

scale = sum(W);
W = W./scale;
H = H.*scale';

options.W = W;
options.H = H;

% TITANized algo
options.inertial = true;
[W0,H0,e,en,ep,etx,lambda] = titanminvol(X,r,options);
titan_error =  norm(X-W0*H0,'fro')/nX;
disp("Titanized relative error = "+titan_error);

%% MM algo

% Run min-vol NMF (default parameters: delta=0.1, maxiter=100)
options = [];
options.lambda = 0.5;
options.maxiter = 1e3;
epsilon = 0.1;
errorX1=[];

for iter = 1:100
    options.iteration = 2;
    lambda1(iter) = options.lambda;
    [W,H] = minvolNMF2(X,r,options);
    options.lambda = 2*sqrt(sum(sum((X - W*H).^2)) + epsilon)*options.lambda;
    options.W = W;
    options.H = H;
    
    errorX1(iter) = norm(X-W*H,'fro')/nX;

    disp("========Iteration "+iter+"=========");
    disp("ErrorX = "+errorX0(iter));
    disp("Lambda = "+lambda1(iter));
    
    
end

W1 = W;
H1 = H;

%% Plot Errors 
set(0, 'DefaultLineLineWidth', 2);
figure;
plot(log10(errorX1),'o-');
xlabel("iteration");
ylabel("log10(errors)");
legend("ErrorX");
title("MM ErrorX");

%% Plot Lambda

set(0, 'DefaultLineLineWidth', 2);
figure;
plot(log10(lambda1),'*-');
xlabel("iteration");
ylabel("log10(Lambdas)");
title("Lambda/iteration");

%% ===2D plot using PCA===
[U,S,V] = svd(X,'econ');

U1 = U(1:2,:);

X_p = U1*X;
W0_p = U1*W0;
W1_p = U1*W1;

set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultLineLineWidth', 2);
figure; 
scatter(X_p(1,:), X_p(2,:),50,'bo'); hold on; 
scatter(W0_p(1,:), W0_p(2,:),100, 'ro','LineWidth', 2); 
scatter(W1_p(1,:), W1_p(2,:),100,'mx','LineWidth', 2); 
grid on; 
legend('Data points $X$','$W^{\ast}$ columns','$\hat{W}$ columns','Interpreter','latex'); 
title("2D plot ($\lambda = 1$)",'Interpreter','latex');




