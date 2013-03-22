function K = MyKernel(X)
% Computes gaussion kernel of the matrix of pairwise distance of the
% samples in / rows of X. sigma is the standard deviation.

n = size(X,1);

D = zeros(n,n);
for i = 1:n
   D(:,i) = sum( (ones(n,1)*X(i,:)-X).^2, 2);
end

sigma = 0.2*max(D(:));
K = exp(- D / (2 * sigma^2));

