% This is the test of the impact of INITIAL LAMBDA in the MM_minvol algo
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
% Goal: Find the largest INITIAL LAMBDA for the MM algorithm that can make
% errorW, and errorX go to the smallest points.
%
% Here, we test the code with each noise level: 0, 0.1, 0.01, 0.001, 0.0001

%% clear data
clear all; clc; 
%close all; 

rng('default')
s = rng

%% synthetic data

% True generating W, with rank_+(W) = 4 > rank(W) = 3
%disp('True W:') 
Wt = [1 0 0 1; 1 0 1 0; 0 1 1 0; 0 1 0 1]' 
m = 15;
r = 6; % nonnegative rank 
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
noiselevel = 0.1;
Xt = Wt*Ht; 
X = max(0,Xt+noiselevel*randn(size(Xt)));
normX = norm(Xt,'fro');

%% Code
lambdaArray = [2,1.5,1,0.8,0.5,0.1,0.05,0.01,0.005,0.001,0.0005,0.0001];
errorWmatrix = [];
errorXmatrix = [];
MMobjmatrix = [];

for i = 1:length(lambdaArray)
% Run min-vol NMF (default parameters: delta=0.1, maxiter=100)

options = [];
errorX0=[];
errorW0=[];
MMobj=[];
options.iteration = 2;
options.lambda = lambdaArray(i);
options.maxiter = 1e3;
epsilon = 0.1;

for iter = 1:100
    lambda1(iter) = options.lambda;
    [W,H,delta] = minvolNMF2(X,r,options);
    
    errorX0(iter) = norm(Xt-W*H,'fro')/normX;
    errorW0(iter) = compareWs(Wt,W);
    MMobj(iter) = sqrt(sum(sum((X - W*H).^2))^2 + epsilon)+options.lambda*log(det(W'*W+delta*eye(size(W'*W))));
    
    %Update lambda,W0,H0
    options.lambda = 2*sqrt(sum(sum((X - W*H).^2)) + epsilon)*options.lambda;
    options.W = W;
    options.H = H;

    disp("========Iteration "+iter+" of lambda = "+lambdaArray(i)+"=========");
    disp("ErrorW = "+errorW0(iter));
    disp("ErrorX = "+errorX0(iter));
    disp("ObjFun = "+MMobj(iter));
    disp("Lambda = "+lambda1(iter));
       
end
errorWmatrix = [errorWmatrix;errorW0];
errorXmatrix = [errorXmatrix;errorX0];
MMobjmatrix = [MMobjmatrix;MMobj];

disp("******* END OF LAMBDA = "+lambdaArray(i)+"**********");
end;

%% Plot ErrorW
set(0, 'DefaultLineLineWidth', 2);
figure;
plot(log10(errorWmatrix(1,:)),'o-');
hold on;
plot(log10(errorWmatrix(2,:)),'+-');
plot(log10(errorWmatrix(3,:)),'*-');
plot(log10(errorWmatrix(4,:)),'.-');
plot(log10(errorWmatrix(5,:)),'x-');
plot(log10(errorWmatrix(6,:)),'s-');
plot(log10(errorWmatrix(7,:)),'d-');
plot(log10(errorWmatrix(8,:)),'^-');
plot(log10(errorWmatrix(9,:)),'v-');
plot(log10(errorWmatrix(10,:)),'>-');
plot(log10(errorWmatrix(11,:)),'<-');
plot(log10(errorWmatrix(12,:)),'p-');
xlabel("iteration",'Interpreter','latex');
ylabel('$\log_{10}$(\textbf{ErrorW})','Interpreter','latex');
grid on;
legend("$\lambda = 2$","$\lambda = 1.5$","$\lambda = 1$","$\lambda = 0.8$","$\lambda = 0.5$","$\lambda = 0.1$","$\lambda = 0.05$","$\lambda = 0.01$","$\lambda = 0.005$","$\lambda = 0.001$","$\lambda = 0.0005$","$\lambda = 0.0001$",'Interpreter','latex');
title("\textbf{ErrorW}",'Interpreter','latex');

%% Plot ErrorX
set(0, 'DefaultLineLineWidth', 2);
figure;
plot(log10(errorXmatrix(1,:)),'o-');
hold on;
plot(log10(errorXmatrix(2,:)),'+-');
plot(log10(errorXmatrix(3,:)),'*-');
plot(log10(errorXmatrix(4,:)),'.-');
plot(log10(errorXmatrix(5,:)),'x-');
plot(log10(errorXmatrix(6,:)),'s-');
plot(log10(errorXmatrix(7,:)),'d-');
plot(log10(errorXmatrix(8,:)),'^-');
plot(log10(errorXmatrix(9,:)),'v-');
plot(log10(errorXmatrix(10,:)),'>-');
plot(log10(errorXmatrix(11,:)),'<-');
plot(log10(errorXmatrix(12,:)),'p-');
xlabel("iteration",'Interpreter','latex');
ylabel('$\log_{10}$(\textbf{ErrorX})','Interpreter','latex');
grid on;
legend("$\lambda = 2$","$\lambda = 1.5$","$\lambda = 1$","$\lambda = 0.8$","$\lambda = 0.5$","$\lambda = 0.1$","$\lambda = 0.05$","$\lambda = 0.01$","$\lambda = 0.005$","$\lambda = 0.001$","$\lambda = 0.0005$","$\lambda = 0.0001$",'Interpreter','latex');
title("\textbf{ErrorX}",'Interpreter','latex');
%% Plot MMObj
set(0, 'DefaultLineLineWidth', 2);
figure;
plot(MMobjmatrix(6,:),'o-');
hold on;
%plot(MMobjmatrix(2,:),'+-');
%plot(MMobjmatrix(3,:),'*-');
%plot(MMobjmatrix(4,:),'.-');
%plot(MMobjmatrix(5,:),'x-');
%plot(MMobjmatrix(6,:),'s-');
plot(MMobjmatrix(7,:),'d-');
plot(MMobjmatrix(8,:),'^-');
plot(MMobjmatrix(9,:),'v-');
plot(MMobjmatrix(10,:),'>-');
plot(MMobjmatrix(11,:),'<-');
plot(MMobjmatrix(12,:),'p-');
xlabel("iteration");
ylabel("MMobj");
legend("0.1","0.05","0.01","0.005","0.001","0.0005","0.0001");
title("MMobj");

%% Plot Lambda

% set(0, 'DefaultLineLineWidth', 2);
% figure;
% plot(laog10(lambda1),'*-');
% xlabel("iteration");
% ylabel("log10(Lambdas)");
% title("Lambda/iteration");