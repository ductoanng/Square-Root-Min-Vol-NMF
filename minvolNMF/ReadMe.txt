The main function, minvolNMF.m, allows you to solve approximately 

min_{W,H} ||M-WH||_F^2 + lambda' * logdet( W^TW + delta I) 

where W >= 0, H >= 0, ||H(:,j)||_1 <= 1 for all j. 

See the paper Minimum-Volume Rank-Deficient Nonnegative Matrix Factorizations, % Valentin Leplat, Andersen M.S. Ang, Nicolas Gillis, 2018. 

You can run synthetic_data to run the same experiment as in that paper. 

