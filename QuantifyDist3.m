function RatioMtrx=QuantifyDist3(DisMtrx1,DisMtrx2,DisMtrx3)

%
results=zeros(3,3);
results(1,3)=1; results(2,3)=1; results(3,3)=1;

% number of realizations
N=size(DisMtrx1,2)-1;
Pyramid = size(DisMtrx1,1);

JS_dispat=DisMtrx1;

JS_CCSIM=DisMtrx2;

JS_sisim=DisMtrx3;



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

if Pyramid ==10
    weight=[1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024]';
else if Pyramid ==3
        weight=[1/2+1/4+1/8+1/16+1/32,1/64+1/128+1/256,1/512+1/1024+1/1024]';
    end
end
% weight=[1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024]';

%weight=[1/2,1/4,1/4]';
ratio_total_disCC=sum(ratio_dis_CC.*weight);
ratio_total_siiCC=sum(ratio_sii_CC.*weight);
ratio_total_dissis=sum(ratio_dis_sii.*weight);
ratio_total_CCsis=sum(ratio_CC_sii.*weight);

results(2,2)=ratio_total_disCC;
results(2,1)=ratio_total_siiCC;

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

if Pyramid ==10
    weight=[1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024]';
else if Pyramid ==3
        weight=[1/2+1/4+1/8+1/16+1/32,1/64+1/128+1/256,1/512+1/1024+1/1024]';
    end
end
% weight=[1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024]';
%weight=[1/2+1/4+1/8+1/16+1/32,1/64+1/128+1/256,1/512+1/1024+1/1024]';
%weight=[1/2,1/4,1/4]';
%weight=[1/2+1/4+1/8+1/16,1/32+1/64+1/128,1/256+1/512+1/1024+1/1024]';
ratio_disCC=sum(ratio_d_C'.*weight);
ratio_siiCC=sum(ratio_s_C'.*weight);
ratio_dissii=sum(ratio_d_s'.*weight);
ratio_CCsii=sum(ratio_C_s'.*weight);
%%

results(1,2)=ratio_disCC;
results(1,1)=ratio_siiCC;
results(3,1)=results(1,1)/results(2,1);
results(3,2)=results(1,2)/results(2,2);
RatioMtrx=results;

end