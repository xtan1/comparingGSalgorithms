function RatioMtrx=QuantifyDist2(DisMtrx1,DisMtrx2)




% Results (3,2):
%                                        DISPAT:CCSIM        CCSIM:CCSIM
%Ratio of between-realization distance     
%Ratio of within-realization distance
%Total ratio
%
results=zeros(3,2);
results(1,2)=1; results(2,2)=1; results(3,2)=1;

% number of realizations
N=size(DisMtrx1,2)-1;
Pyramid = size(DisMtrx1,1);

JS_dispat=DisMtrx1;

JS_CCSIM=DisMtrx2;





d_dispat=sum(JS_dispat,3);
d_CCSIM=sum(JS_CCSIM,3);
% sum of distances between training image and all 50 other realizations
d_disPR=d_dispat(:,51);
d_CCSIMPR=d_CCSIM(:,51);


ratio_dis_CC=d_disPR./d_CCSIMPR;
ratio_CC_dis=d_CCSIMPR./d_disPR;


if Pyramid ==10
    weight=[1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024]';
else if Pyramid ==3
        weight=[1/2+1/4+1/8+1/16+1/32,1/64+1/128+1/256,1/512+1/1024+1/1024]';
    end
end
% weight=[1/2,1/4,1/8,1/16,1/32,1/64,1/128,1/256,1/512,1/1024]';

%weight=[1/2,1/4,1/4]';
ratio_total_disCC=sum(ratio_dis_CC.*weight);
ratio_total_CCdis=sum(ratio_CC_dis.*weight);

results(2,1)=ratio_total_disCC;


%% Distance between realizations

JS_dispat50 =JS_dispat(:,1:N,1:N);
JS_CCSIM50 =JS_CCSIM(:,1:N,1:N);



%%

nn=1;
for i=1:Pyramid
    for j=1:50
        for k=1:50
            if j>k
                array1(nn)=JS_dispat50(i,j,k);
                array2(nn)=JS_CCSIM50(i,j,k);
                
                nn=nn+1;
            end 
            
        end
    end
end
%%
dispat10_1225=reshape(array1,1225,Pyramid);
CCSIM10_1225=reshape(array2,1225,Pyramid);


%%

dispat10_1225_avg=mean(dispat10_1225);
CCSIM10_1225_avg=mean(CCSIM10_1225);

ratio_d_C=dispat10_1225_avg./CCSIM10_1225_avg;

ratio_C_d=CCSIM10_1225_avg./dispat10_1225_avg;


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


ratio_CCdis=sum(ratio_C_d'.*weight);
%%

results(1,1)=ratio_disCC;
results(3,1)=results(1,1)/results(2,1);

RatioMtrx=results;
% Results (3,2):
%                                        DISPAT:CCSIM        CCSIM:CCSIM
%Ratio of between-realization distance     
%Ratio of within-realization distance
%Total ratio
%

end