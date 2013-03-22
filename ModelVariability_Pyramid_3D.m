%% step 2
%% 3D case

clc;
clear all;

N=50;
Pyramid = 10; % depend on the size of realizations and the size of patterns

%% MR
load 'Ti_3Dchannel.mat';
% 1. for diapat
%load 're_bi_3DDispat.mat';
% 2. for CCSIM 
load 'CCSIM_3D_binary.mat';
% 3. for SNESIM
%load 're_bi_3Dsnesim.mat';

realization1=reshape(Final,69,69,39,50);
%load 'i_channel.mat';

fp1      = mfilename('fullpath');
dirName1 = fileparts(fp1);
slash1   = strfind(dirName1, '\');
dirName1 = dirName1(1:slash1(end)-1);
addpath([dirName1 '\Copy of SimplifiedProject\dist_pat']);

dim = size(out);

%% Pyramid = 10; realization = 101*101
hist=-1*ones(Pyramid, N+1, 262144);
%prob=-1*ones(Pyramid, N, 65535);
for MR=1:10
    for i=1:N
    out_c=realization1(:,:,:,i);
    out1 = my_resize(out_c,ceil(dim(1)/MR),ceil(dim(2)/MR),ceil(dim(3)/MR)); 
    level = graythresh(out1);
    out1 = im2bw_3D(out1, level);
    temp=calculateMPH_3D(out1);
    hist(MR,i,:)=temp(1:262144);
    hist(MR,i,:)=hist(MR,i,:)/sum(hist(MR,i,:));
    end
    out_c=out;
    out1 = my_resize(out_c,ceil(dim(1)/MR),ceil(dim(2)/MR),ceil(dim(3)/MR)); 
    level = graythresh(out1);
    out1 = im2bw_3D(out1, level);
    temp=calculateMPH_3D(out1);
    hist(MR,N+1,:)=temp(1:262144);
    hist(MR,N+1,:)=hist(MR,N+1,:)/sum(hist(MR,N+1,:));
    
end





%% calculate JS divergence among hist

JS_Pyramid3D=-1*ones(Pyramid,N+1,N+1);
X=1:262144;

for MR=1:Pyramid
    for i=1:N+1
        for j=1:N+1 
            temp1=hist(MR,i,:);temp2=hist(MR,j,:);
            temp1=temp1(:);temp2=temp2(:);
            JS_Pyramid3D(MR,i,j) = kldiv(X,temp1'+eps,temp2'+eps,'js');
            
        end
    end
end


%save('JS_sisim_MR.mat', 'JS1_sisim','JS2_sisim','JS3_sisim');
%save('JS_CCSIMMV3D.mat', 'JS_Pyramid3D');
save('JS_SNESIMMV3D.mat', 'JS_Pyramid3D');
%save('JS_DISPATMV3D.mat', 'JS_Pyramid3D');