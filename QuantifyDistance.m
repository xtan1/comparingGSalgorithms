%% step 4


clc;
clear all;
close all;

% number of realizations
N=50;
Pyramid = 10;
% 2D case
load JS_dispatMV.mat;
JS_dispat=JS_Pyramid;
load JS_CCSIMMV.mat;
JS_CCSIM=JS_Pyramid2;
load JS_sisim_MV.mat;
JS_sisim=JS_Pyramid3;

% SNESIM 101*101 15*15
% load JS_SNESIMMV_1515.mat;
% JS_dispat=JS_Pyramid_15_15;
% load JS_SNESIMMV_3515.mat;
% JS_CCSIM=JS_Pyramid_35_15;
% load JS_SNESIMMV_40015.mat;
% JS_sisim=JS_Pyramid_400_15;

% 2D continuous case 260x200
% load JS_dispatMV_cont.mat;
% JS_dispat=JS_Pyramid2;
% load JS_CCSIM130260MV_cont.mat;
% JS_CCSIM=JS_Pyramid2;
% load 'JS_SGSIM130260MV_cont.mat';
% JS_sisim=JS_Pyramid2;

% load full_cluster_DISPAT_new.mat;
% JS_dispat=JS_Pyramid2;
% load full_cluster_CCSIM_new.mat;
% JS_CCSIM=JS_Pyramid2;
% load full_cluster_SISIM_new.mat;
% JS_sisim=JS_Pyramid2;

% SNESIM 101*101 101*101
% load JS_SNESIMMV_10101.mat;
% JS_dispat=JS_Pyramid_10_101;
% load JS_SNESIMMV_50101.mat;
% JS_CCSIM=JS_Pyramid_50_101;
% load JS_SNESIMMV_400101.mat;
% JS_sisim=JS_Pyramid_400_101;

% 3D case
% load JS_DISPATMV3D.mat;
% JS_dispat=JS_Pyramid3D;
% load JS_CCSIMMV3D.mat;
% JS_CCSIM=JS_Pyramid3D;
% load JS_SNESIMMV3D.mat;
% JS_sisim=JS_Pyramid3D;

% % 2D continuous case
% load JS_dispatMV_cont.mat;
% JS_dispat=JS_Pyramid2;
% load JS_CCSIMMV_cont.mat;
% JS_CCSIM=JS_Pyramid2;
% load JS_SGSIMMV_cont.mat;
% JS_sisim=JS_Pyramid2;

d_dispat=sum(JS_dispat,3);
d_CCSIM=sum(JS_CCSIM,3);
d_sisim=sum(JS_sisim,3);
% sum of distances between training image and all 50 other realizations
d_disPR=d_dispat(:,51);
d_CCSIMPR=d_CCSIM(:,51);
d_sisimPR=d_sisim(:,51);

ratio_dis_CC=d_disPR./d_CCSIMPR;
ratio_sii_CC=d_sisimPR./d_CCSIMPR;
ratio_CC_sii=d_CCSIMPR./d_sisimPR;
ratio_dis_sii=d_disPR./d_sisimPR;
weight=[1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024]';
% weight=[1/2+1/4+1/8+1/16+1/32,1/64+1/128+1/256,1/512+1/1024+1/1024]';
%weight=[1/2,1/4,1/4]';
ratio_total_disCC=sum(ratio_dis_CC.*weight);
ratio_total_siiCC=sum(ratio_sii_CC.*weight);
ratio_total_dissis=sum(ratio_dis_sii.*weight);
ratio_total_CCsis=sum(ratio_CC_sii.*weight);
%% Distance between realizations

JS_dispat50 =JS_dispat(:,1:N,1:N);
JS_CCSIM50 =JS_CCSIM(:,1:N,1:N);
JS_sisim50 =JS_sisim(:,1:N,1:N);


%%

nn=1;
for i=1:Pyramid
    for j=1:50
        for k=1:50
            if j>k
                array1(nn)=JS_dispat50(i,j,k);
                array2(nn)=JS_CCSIM50(i,j,k);
                array3(nn)=JS_sisim50(i,j,k);
                nn=nn+1;
            end 
            
        end
    end
end
%%
dispat10_1225=reshape(array1,1225,Pyramid);
CCSIM10_1225=reshape(array2,1225,Pyramid);
sisim10_1225=reshape(array3,1225,Pyramid);

%%

dispat10_1225_avg=mean(dispat10_1225);
CCSIM10_1225_avg=mean(CCSIM10_1225);
sisim10_1225_avg=mean(sisim10_1225);

ratio_d_C=dispat10_1225_avg./CCSIM10_1225_avg;
ratio_s_C=sisim10_1225_avg./CCSIM10_1225_avg;
ratio_C_s=CCSIM10_1225_avg./sisim10_1225_avg;
ratio_d_s=dispat10_1225_avg./sisim10_1225_avg;

weight=[1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024]';
% weight=[1/2+1/4+1/8+1/16+1/32,1/64+1/128+1/256,1/512+1/1024+1/1024]';
%weight=[1/2,1/4,1/4]';
%weight=[1/2+1/4+1/8+1/16,1/32+1/64+1/128,1/256+1/512+1/1024+1/1024]';
ratio_disCC=sum(ratio_d_C'.*weight);
ratio_siiCC=sum(ratio_s_C'.*weight);
ratio_dissii=sum(ratio_d_s'.*weight);
ratio_CCsii=sum(ratio_C_s'.*weight);
%%


%% Distance between ti and realization

JS_dispat51=squeeze(JS_dispat(:,N+1,1:N));
JS_CCSIM51 =squeeze(JS_CCSIM(:,N+1,1:N));
JS_sisim51 =squeeze(JS_sisim(:,N+1,1:N));

dispat51_sum=sum(JS_dispat51,2);
CCSIM51_sum=sum(JS_CCSIM51,2);
sisim51_sum=sum(JS_sisim51,2);

dispat51_avg=dispat51_sum/50;
CCSIM51_avg=CCSIM51_sum/50;
sisim51_avg=sisim51_sum/50;

total_dispat51_avg=sum(dispat51_avg.*weight);
total_CCSIM51_avg=sum(CCSIM51_avg.*weight);
total_sisim51_avg=sum(sisim51_avg.*weight);

% use variance
dispat51_std=std(JS_dispat51');
CCSIM51_std=std(JS_CCSIM51');
sisim51_std=std(JS_sisim51');

total_dispat51=sum(dispat51_std'.*weight);
total_CCSIM51=sum(CCSIM51_std'.*weight);
total_sisim51=sum(sisim51_std'.*weight);
