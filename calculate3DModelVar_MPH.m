function DistMtrx = calculate3DModelVar_MPH(realization1,out,pyramid)
    % This function is suitable to 3D binary case
    % input: realization1 is all the realizations needed: e.g.69*69*39*50
    % input: out is the training image: e.g. 69*69*39
    % output: DistMtrx: e.g. 10*51*51 (10 is the level of Pyramid)

    N=size(realization1,4); % number of realizations
    Pyramid = pyramid; % depend on the size of realizations and the size of patterns



    dim = size(out);

    %% Pyramid = 10; realization = 101*101
    % 262144=2^(3*3*2): 3*3*2 is the size of the template in 3D
    hist=-1*ones(Pyramid, N+1, 262144,'single');
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
   
    hist3=hist;
    DistMtrx=JS_Pyramid3D;

end