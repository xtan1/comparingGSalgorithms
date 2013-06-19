% test
clc;
clear all;
close all;

fp1      = mfilename('fullpath');
dirName1 = fileparts(fp1);
slash1   = strfind(dirName1, '\');
dirName1 = dirName1(1:slash1(end)-1);
addpath([dirName1 '\ComparingGSAlgorithms\Data']);

%% 2D binary case
% load 'SISIM_channel.mat';
% load 'Ti_channel.mat';
% % the first parameter is the realizations
% % the second parameter is the training image
% % the third parameter is the number of resolution
% DistanceMatrix1=calculateModelVar_MPH(realization1,out,10);

%% 3D binary case
clear all;
load 'Ti_3Dchannel.mat';
load 'CCSIM_3D_binary.mat';
realization1=reshape(Final,69,69,39,50);
% % easy to be out of memory, may need to preallocate the memory
% DistanceMatrix2=calculate3DModelVar_MPH(realization1,out,10);

%% continuous case
% clear all;
% load 'TI_cont_ISO.mat';
% out=realization1;
% 
% %load 'new_cont.mat';
% load 'CCSIM_cont_re.mat';
% %realization1=new_ccsim;
% tempSize =23;
% DistMtrx2 = calculateModelVar_CHP(realization1,out,tempSize);
% save('JS_CCSIM_lasttry.mat', 'DistMtrx2');
% 
% 
% clear all;
% load 'TI_cont_ISO.mat';
% out=realization1;
% 
% load 'SGS_cont_re.mat';
% tempSize =23;
% DistMtrx1 = calculateModelVar_CPH(realization1,out,tempSize);
% save('JS_SGSIMlasttry.mat', 'DistMtrx1');


%% 3D case
clear all;
load 'Ti_3Dchannel.mat';
load 'CCSIM_3D_binary.mat';
realization1=reshape(Final,69,69,39,50);

% 
% %load 'new_cont.mat';
% load 'CCSIM_cont_re.mat';
tempSize =[23,23,10];
% the number of patterns you will skip when scanning TI
% increase if out of memory
skipPattern=[8,8,6];
DistMtrx2 = calculate3DModelVar_CHP(realization1,out,tempSize,skipPattern);
save('JS_CCSIM_finalfinal.mat', 'DistMtrx2');




%% show 2D case 
% clear all;
% load full_cluster_DISPAT_new.mat;
% DisMtrx1=JS_Pyramid2;
% load full_cluster_CCSIM_new.mat;
% DisMtrx2=JS_Pyramid2;
% load full_cluster_SISIM_new.mat;
% DisMtrx3=JS_Pyramid2;
% 
% load 'Ti_channel.mat';
% % realizations
% load 'realization_dispat.mat';
% re1=realization1;
% load 'CCSIM_new.mat';
% re2=reshape(Final,101,101,50);
% load 'sisim_channel.mat';
% re3=realization1;

%show2DModelVar2(DisMtrx1,DisMtrx2,out,re1,re2);
% show2DModelVar3(DisMtrx1,DisMtrx2,DisMtrx3,out,re1,re2,re3);
%% show 3D case
% clear all;
% load JS_DISPATMV3D.mat;
% JS_dispat=JS_Pyramid3D;
% load JS_CCSIMMV3D.mat;
% JS_CCSIM=JS_Pyramid3D;
% load JS_SNESIMMV3D.mat;
% JS_sisim=JS_Pyramid3D;
% 
% load 'Ti_3Dchannel.mat';
% load 're_bi_3DDispat.mat'
% re1=realization1;
% load 're_bi_3Dsnesim.mat';
% re3=realization1;
% load 'CCSIM_3D_binary.mat'
% model_CCSIM=reshape(Final,69,69,39,50);
% re2=model_CCSIM;
% 
%show3DModelVar2(JS_dispat,JS_CCSIM,out,re1,re2);
% show3DModelVar3(JS_dispat,JS_CCSIM,JS_sisim,out,re1,re2,re3);
%% Quantify Distance
% clear all;
% % 
% % % 2D case
% load JS_dispatMV.mat;
% JS_dispat=JS_Pyramid;
% load JS_CCSIMMV.mat;
% JS_CCSIM=JS_Pyramid2;
% load JS_sisim_MV.mat;
% JS_sisim=JS_Pyramid3;

%RatioMtrx2=QuantifyDist2(JS_dispat,JS_CCSIM);
% RatioMtrx3=QuantifyDist3(JS_dispat,JS_CCSIM,JS_sisim);

% 3D case
% clear all;
% load JS_DISPATMV3D.mat;
% JS_dispat=JS_Pyramid3D;
% load JS_CCSIMMV3D.mat;
% JS_CCSIM=JS_Pyramid3D;
% RatioMtrx2=QuantifyDist2(JS_dispat,JS_CCSIM);
