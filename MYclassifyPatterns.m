function [X, Y, e, K, idx, prototype, MDS] = MYclassifyPatterns(rea, par)


% ------------ INITIALIZATION ------------
Patx = par.Patx;
Paty = par.Paty;
Patz = par.Patz;
sN   = Patx*Paty*Patz;
nr   = size(rea,4);
sV   = max(rea(:));

Dimx = par.Dimx;
Dimy = par.Dimy;
Dimz = par.Dimz;
mx   = par.mx;
my   = par.my;
mz   = par.mz;

m1   = par.m1;
clus = par.clus;



% ------------ EXTRACT PATTERNS ------------ 
disDimx = ceil((Dimx - (1+(Patx-1)*m1) + 1)/mx);
disDimy = ceil((Dimy - (1+(Paty-1)*m1) + 1)/my);
%%

%disDimx = Dimx-Patx +1;
%disDimy = Dimy-Paty +1;
%%
if Dimz > 1
    disDimz = ceil((Dimz - (1+(Patz-1)*m1) + 1)/mz);
    %disDimz = Dimz-Patz +1;
else
    disDimz = 1;
end
   
l=0;
for i=1:disDimx
    for j=1:disDimy
        for k=1:disDimz
            
            wx = 1+mx*(i-1):m1:1+mx*(i-1)+(Patx-1)*m1;
            %wx = i:m1:i+(Patx-1);
            wy = 1+my*(j-1):m1:1+my*(j-1)+(Paty-1)*m1;
            %wy = j:m1:j+(Paty-1);
            if Dimz > 1
                wz = 1+mz*(k-1):m1:1+mz*(k-1)+(Patz-1)*m1;
                %wz = k:m1:k+(Patz-1);
            else
                wz = 1;
            end

            for s=1:nr
                l=l+1;
                X(l,:)=reshape(rea(wx,wy,wz,s),1,Patx*Paty*Patz);               
            end
        end
    end
end
fprintf('Patterns retained for analysis %d = %d x %d.\n', l, nr, l/nr);



% ------------ DISSIMILARITY MATRIX ------------ 
fprintf('Dissimilarity Matrix is being built');

% If binary image 
if all(X(:) == 0 | X(:) == 1)
   way = 2;                      % Use distance transform (proximity) 
   fprintf(' using the distance transform... ')
else
   way = 1;                      % Use euclidean distance
   fprintf(' using the Euclidean distance... ')
end
disMat = compDisMat(X, Patx, Paty, Patz, way);
fprintf(' Done!\n');

if 0
    figure
    imagesc(disMat)
    axis image
    colorbar
end

% ------------ MULTI DIMENSIONAL SCALING ------------
% Select intrinsic dimensionality using PCA

if isfield(par, 'MDS')
   MDS = par.MDS;
else
   fprintf('Calculating Intrinsic Dimensionality ..........');
   S   = slpca(X','preserve',{'Everrit',0});   % other method can be profileML
   MDS = S.feadim;
   fprintf('\b\b\b\b %3d\n',MDS)
end

fprintf('MultiDimensional Scaling of dimension = %d.', MDS);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
   [Y1, e] = cmdscale(disMat);
catch
   warning('Uses MYcmdscale instead of cmdscale')
   [Y1, e] = MYcmdscale(disMat);
   
end
Y = Y1(:,1:min(MDS,end));

fprintf(' Done!\n');


% ------------ CLUSTERING ALGORITHM ------------ 
fprintf('Kernel K-means analysis.');

% Apply Radial Basis Function Kernel
K = MyKernel(Y);

% Kernel K-means algorithm
idx = dualkmeansFast(K, clus);

% Computing prototypes
prototype = zeros(clus, sN);
j = 0; 
for i = 1:clus
   clust = find(idx==i);
   n = length(clust);
   
   if n > 0
      j = j + 1;
      idx(clust) = j;
      
      prototype(j,:)=sum(X(clust,:),1)./n;
      
      % Display prototypes and cluster members
      if 0 %Patz == 1
         subplot(1,2,1)
         imagesc(reshape(prototype(j,:), Patx, Paty))
         axis image
         colormap gray
         caxis([0 sV])
         
         for k = 1:n
            subplot(1,2,2)
            imagesc(reshape( X(clust(k),:), Patx, Paty))
            axis image
            pause(0.5)
            caxis([0 sV])
         end
      end
   end
end
prototype = prototype(1:j,:);

fprintf(' Done!\n');



function D = compDisMat(X, Patx, Paty, Patz, way)

npat = size(X,1);

switch way
   
   % Euclidean
   case 1
      D = zeros(npat, npat);
      
      % Use euclidean distance
      for i = 1 : npat
         for j = i+1 : npat
            
            dX = X(i,:) - X(j,:);
            D(i,j) = sqrt(dX * dX');
            D(j,i) = D(i,j);
         end
      end
      
   % Distance transform (proximity)   
   case 2
      
      sN = size(X,2);
      temp=zeros(npat,sN);
      for i = 1 : npat
         DT=bwdist(reshape(X(i,:), Patx, Paty, Patz));
         if DT==Inf
            DT = ones(1,sN);
         end
         
         % KATLA: changed as all foreground patterns became NaN
         maxDT = max(DT(:));
         if maxDT > 0
             temp(i,:)=reshape(1-DT./maxDT,1,[]);
         else
             temp(i,:)=1;
         end
      end
      
      D = compDisMat(temp, Patx, Paty, Patz, 1);
end