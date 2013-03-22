function [H, n] = CompHistRea(rea, prototype, par, mindist)

% Dimension of patterns
nZ = par.Dimx; mZ = par.Dimy; pZ = par.Dimz;
nN = par.Patx; mN = par.Paty; pN = par.Patz;
dn = par.mx; dm = par.my; dp = par.mz;
m1 = par.m1;
ncat = size(prototype,1);
E = ones(ncat,1);

% Histogram
H = zeros(ncat,1);

% No minimum distance
if nargin < 4
   mindist = inf;
end

% For all inner voxels
n = 0;
for i = 1 : dn : nZ-m1*(nN-1)
   dim1 = i : m1 : i+m1*(nN-1);
   
   for j = 1 : dm : mZ-m1*(mN-1)
      dim2 = j : m1 : j+m1*(mN-1);
         
      for k = 1 : dp : pZ-m1*(pN-1)
         dim3 = k : m1 : k+m1*(pN-1);

         n = n + 1;
         
         % Voxel values of the pattern
         zijk = rea(dim1,dim2,dim3);
         
         % Find distances to all prototypes
         dist = sum(abs(prototype - E*zijk(:)'),2);
         
         % Find minimum distance and thereby pattern category
         [~, id] = min(dist);
         
         % Update the frequency distribution if sufficiently close
         if dist(id) <= mindist
            H(id) = H(id) + 1;
         end
      end
   end
end
