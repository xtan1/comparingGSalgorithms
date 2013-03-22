function d = ChiDist(H1, H2, n1, n2)

ids = H1 > 0 | H2 > 0;

H1 = H1/n1;
H2 = H2/n2;
n1 = 1;
n2 = 1;

% Compute the distance
d = 1./(H1(ids)+H2(ids))' * (sqrt(n1/n2)*H2(ids) - sqrt(n2/n1)*H1(ids)).^2 ...
    + n1 - n1/n2 * sum(H2);  
   
 if isnan(d)
    keyboard
 end