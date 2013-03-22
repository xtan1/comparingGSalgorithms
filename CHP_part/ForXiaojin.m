clear all
close all
clc

% Load realizations
%load saturation_profile3.mat
load 'CCSIM_2D_Cont.mat';
sat_prof = reshape(C,130,100,50);
nr = size(sat_prof,3);

parS.clus = 50;

% Visualize data
nxS = size(sat_prof,1);
nyS = size(sat_prof,2);
nzS = 1;
%%
% if 1
%    figure
%    for i = 1:nr
%       imagesc(sat_prof(:,:,i))
%       colorbar
%       caxis([0 1])
%       axis image
%       title(sprintf('Realization %d out of %d', i, nr))
%       
%       pause(0.1)
%    end
% end
%%

% Dimensions
parS.multipleGrid = 1;
parS.m1 = 1;
parS.MDS = 20;

reaS = reshape(sat_prof,[nxS,nyS,nzS,nr]);
parS.Patx = 23;
parS.Paty = 23;
parS.Patz = 1;
%parS.mx = parS.Patx;
%parS.my = parS.Paty;
%parS.mz = parS.Patz;
parS.mx=2;
parS.my=2;
parS.mz=1;
parS.Dimx = nxS;
parS.Dimy = nyS;
parS.Dimz = nzS;

% Training image
%TIS = reaS;
TIS=imin;
% We dont have a TI so we use the realizations themselves.

% Classify patterns
[XS, YS, eS, ~, idxS, prototypeS, MDSS] = MYclassifyPatterns(TIS, parS);
parS.clus = size(prototypeS,1);


% Build the clustered histograms
HS = zeros(parS.clus,nr);

parS.mx = 1;
parS.my = 1;
parS.mz = 1;

for i = 1:nr
   [HS(:,i), npatS] = CompHistRea(reaS(:,:,i), prototypeS, parS);
end

% Compute distances between clustered histograms of models
DS = zeros(nr, nr);

for i = 1:nr
   for j = i+1:nr     
      %DS(i,j) = ChiDist(HS(:,i), HS(:,j), npatS, npatS);
      DS(i,j) = JSDist(HS(:,i), HS(:,j));
      DS(j,i) = DS(i,j);
   end
end


% Plots
figure
imagesc(DS)
title('Distance matrix for the saturation models')
axis image
colorbar

