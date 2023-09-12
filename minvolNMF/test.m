% % X = readmatrix("CreditCardCustomerData.csv");
% % X(:,1) = [];
% % %disp(X);readmatrix("CreditCardCustomerData.csv")
% % X = normalize(X,"norm"); 
% % disp([X X]);
% 
% % M = [2 3 4 1 4];
% % 
% % figure;
% % plot(M,'ro--')
% 
% 
% %Wt = [1 0 0 1; 1 0 1 0; 0 1 1 0; 0 1 0 1]';
% 
% %test
% Wt = [1 3 2 4; 0 1 5 3; 0 0 2 2; 0 0 0 1]
% 
% Wt(1:2,:)
% [Wt;Wt(1:2,:)]
% 
% r = 4; % nonnegative rank 
% % Generate H with n columns
% n = 500; 
% purity = 0.8; % ~min distance between the columns of W and the data points
% alpha = 0.05*ones(r,1);  % Parameter of the Dirichlet distribution
% Ht = [sample_dirichlet(alpha,n)']; 
% for j = 1 : n
%     while max( Ht(:,j) ) > purity
%         Ht(:,j) = sample_dirichlet(alpha,1)';
%     end
% end
% 
% % Xt = Wt*Ht;
% % W1 = mvsa(Xt,r);
% % disp("Simplex error: "+compareWs(Wt,W1));
% 
% listK = [1.5 1 0.5 0.1 0.09 0.05 0.01 0.005 0.001 0.0009 0.0007 0.0005 0.0004 0.00035 0.0003 0.00025 0.0002 0.0001 0.00005 0.00001 0.000005 0.000001 1e-07 1e-08 1e-09 1e-10 1e-12 1e-15];
% for k = 1:length(listK)
%     options.lambda=listK(k);
%     %disp("lambda~ = "+ options.lambda);
%     [W,H,e,er1,er2,lambda] = minvolNMF(Xt,r,options);
%     %disp("lambda = "+lambda)
%     trueLambda(k) = options.lambda;
%     %Wmatrix(k) = W;
%     %Error ||W-Wt||/||Wt||
%     errorW(k) = compareWs( Wt, W );
%     errorX(k) = compareWs( X,W*H);
%     %fprintf('Error ||W-Wt||/||Wt|| = %2.2f%%.\n', error(k) ); 
% end
% 
% [We_noi,indexW] = min(errorW);
% options.lambda=listK(indexW);
% [W,H,e,er1,er2,lambda] = minvolNMF(Xt,r,options);
% %[Xe_noi(j),indexX] = min(errorX);
% disp("PFGD Min Vol ||W-Wt|| error : "+ We_noi);
% disp("Compare W and W': "+compareWs(W1,W));
% 
% 
clear all;
X = randperm(4,4)
X = [0 0 0 1; 1 0 0 0; 0 0 1 0; 0 1 0 0]

Y = [1 2 3 4; 2 4 6 8; 7 8 9 10; 11 12 13 14]
rank(Y)
Z=X*Y
rank([Y;Z])
