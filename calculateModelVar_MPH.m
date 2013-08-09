function DistMtrx = calculateModelVar_MPH(realization1,out,pyramid)
    % input: realization1 is all the realizations needed: e.g.101*101*50
    % input: out is the training image
    % inout: pyramid is the level of Pyramid
    % output: DistMtrx: e.g. 10*51*51 (10 is the level of Pyramid)

    N=size(realization1,3); % number of realizations
    Pyramid = pyramid; % depend on the size of realizations and the size of patterns


    %% Pyramid = 10; realization = 101*101
    % 65536=2^(4*4): 4*4 is the size of the template
    hist=-1*ones(Pyramid, N+1, 65536);
    %prob=-1*ones(Pyramid, N, 65535);
    for MR=1:10
        for i=1:N
        out_c=realization1(:,:,i);
        out1 = imresize(out_c,1/MR); 
        level = graythresh(out1);
        out1 = im2bw(out1, level);
        temp=calculateMPH(out1);
        hist(MR,i,:)=temp(1:65536);
        hist(MR,i,:)=hist(MR,i,:)/sum(hist(MR,i,:));
        end
        out_c=out;
        out1 = imresize(out_c,1/MR); 
        level = graythresh(out1);
        out1 = im2bw(out1, level);
        temp=calculateMPH(out1);
        hist(MR,N+1,:)=temp(1:65536);
        hist(MR,N+1,:)=hist(MR,N+1,:)/sum(hist(MR,N+1,:));
    
    end





    %% calculate JS divergence among hist

    JS_Pyramid=-1*ones(Pyramid,N+1,N+1);
    X=1:65536;

    for MR=1:Pyramid
        for i=1:N+1
            for j=1:N+1 
                temp1=hist(MR,i,:);temp2=hist(MR,j,:);
                temp1=temp1(:);temp2=temp2(:);
                JS_Pyramid(MR,i,j) = kldiv(X,temp1'+eps,temp2'+eps,'js');
            
            end
        end
    end

    JS_Pyramid3=JS_Pyramid;
    %save('JS_sisim_MR.mat', 'JS1_sisim','JS2_sisim','JS3_sisim');
    hist3=hist;
    DistMtrx=JS_Pyramid3;

end
