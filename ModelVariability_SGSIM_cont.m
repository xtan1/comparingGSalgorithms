clear all
close all
clc

fp1      = mfilename('fullpath');
dirName1 = fileparts(fp1);
slash1   = strfind(dirName1, '\');
dirName = dirName1(1:slash1(end));
dirName1 = dirName1(1:slash1(end)-1);
addpath([dirName1 '\Copy of SimplifiedProject\dist_pat']);
addpath([dirName 'New Code\CHP_part']);
% Load realizations
%load saturation_profile3.mat

% obtain training image
% load 'CCSIM_2D_Cont12.mat';
% clear C;
%% 260*200 case
load 'sgsim_260_200.mat';
imin=realization1;
%%


% %load 're_cont_2DDispat.mat';
% load 're_cont_2DSGSIM_original.mat';
% sat_prof = reshape(realization1,130,100,50);
%% 260*200 case
load 're_sgsim_260_200.mat';
sat_prof = reshape(realization1,260,200,50);
%%


nr = size(sat_prof,3);

parS.clus = 50;

Pyramid = 3;

DS = zeros(Pyramid, nr+1, nr+1);


%%
% for iii=1:50
%     figure;
%     imshow(sat_prof(:,:,iii));
%     colormap('jet');
% end
% %

HHH_SGSIM=zeros(Pyramid,parS.clus,nr+1);
% pts=zeros(Pyramid,parS.clus,23*23);

for MR=1:Pyramid
    
    
    
    TIS1=imresize(imin,1/MR);
    
    % Visualize data
    nxS = size(TIS1,1);
    nyS = size(TIS1,2);
    nzS = 1;
    
    % Dimensions
    parS.multipleGrid = 1;
    parS.m1 = 1;
    parS.MDS = 20;
    
    for ii=1:nr
        aaa(:,:,ii)=imresize(sat_prof(:,:,ii),1/MR);
    end

    %reaS = reshape(sat_prof,[nxS,nyS,nzS,nr]);
    reaS = reshape(aaa,[nxS,nyS,nzS,nr]);
    
    %
    clear aaa;
    %
    parS.Patx = 43;
    parS.Paty = 43;
    parS.Patz = 1;
    %parS.mx = parS.Patx;
    %parS.my = parS.Paty;
    %parS.mz = parS.Patz;
    parS.mx=4;
    parS.my=4;
    parS.mz=1;
    parS.Dimx = nxS;
    parS.Dimy = nyS;
    parS.Dimz = nzS;
    
    
    % Classify patterns
    [XS, YS, eS, ~, idxS, prototypeS, MDSS] = MYclassifyPatterns(TIS1, parS);
    parS.clus = size(prototypeS,1);
    
    % Bulid the histogram of training image
    hh=zeros(parS.clus,1);
    for iii = 1:size(idxS,1)
        hh(idxS(iii))=hh(idxS(iii))+1;
    end
    
    % Dispaly the prototypes and its counts
%     figure;
%     for mm=1:parS.clus
%         
%         subplot(7,7,mm)
%         pt=reshape(prototypeS(mm,:),[parS.Patx,parS.Paty]);
%         imshow(pt);
%         colormap('jet');
% %         title(['Ti, MR=' int2str(MR), ' Count=', int2str(hh(mm))]);
%         title(int2str(hh(mm)));
%     end
    
    % Build the clustered histograms
    HS = zeros(parS.clus,nr+1);

    parS.mx = 1;
    parS.my = 1;
    parS.mz = 1;

    for i = 1:nr
        [HS(:,i), npatS] = CompHistRea(reaS(:,:,i), prototypeS, parS);
%         figure;
%         for mm=1:parS.clus
%         
%             subplot(7,7,mm)
%             pt=reshape(prototypeS(mm,:),[parS.Patx,parS.Paty]);
%             imshow(pt);
%             colormap('jet');
%             title(['re', int2str(i), 'MR=' int2str(MR), ' Count=', int2str(HS(mm,i))]);
%         end
    end
    HS(:,nr+1) = hh;
    NClust=size(HS,1);
    HHH_SGSIM(MR,1:NClust,:)=HS;
    % Store prototypeS
%     pts(MR,1:NClust,:)=prototypeS;
    
    % Compute distances between clustered histograms of models
    
    %for MRT=1:Pyramid
        for i = 1:nr+1
            for j = i+1:nr+1     
                %DS(i,j) = ChiDist(HS(:,i), HS(:,j), npatS, npatS);
                DS(MR,i,j) = JSDist(HS(:,i), HS(:,j));
                DS(MR,j,i) = DS(MR,i,j);
            end
        end 
    %end
    
    
end


JS_Pyramid2=DS;
%save('JS_sisim_MR.mat', 'JS1_sisim','JS2_sisim','JS3_sisim');
save('JS_SGSIM260MV_cont.mat', 'JS_Pyramid2','HHH_SGSIM');


