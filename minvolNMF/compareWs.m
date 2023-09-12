% Compute
%
% errW = min_perm ||W - West(:,perm)||_F / ||W||_F

function errW = compareWs( W, West );

r = size(W,2); 

for i = 1 : r
    for j = 1 : r
        Dist(i,j) = norm( W(:,i) - West(:,j) , 'fro' )^2; 
    end
end
 
[assignment,cost] = munkres(Dist);
%assignment is an array of column after permutation 
%by Munkres Assignment Algorithm for lowest cost

errW = norm( W - West(:,assignment) , 'fro' ) / norm(W,'fro'); 