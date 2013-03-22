%% step 3
% 3D case


clc;
clear all;
close all;
%% with training image
%% MDS
load JS_DISPATMV3D.mat;
JS_dispat=JS_Pyramid3D;
load JS_CCSIMMV3D.mat;
JS_CCSIM=JS_Pyramid3D;
load JS_SNESIMMV3D.mat;
JS_sisim=JS_Pyramid3D;

%% training image
%load 'Ti_channel.mat';
load 'Ti_3Dchannel.mat';
%% realizations
%load 'realization_dispat.mat';
load 're_bi_3DDispat.mat'
model_dispat=realization1;
model_dispat(:,:,:,51)=out;


load 're_bi_3Dsnesim.mat';
model_sisim=realization1;
model_sisim(:,:,:,51)=out;

%for i=1:50
%    realization1(:,:,i)=realization1(:,:,i)';
%end
%model_sisim=realization1;
%model_sisim(:,:,51)=out;

%load 'CCSIM_new.mat';
load 'CCSIM_3D_binary.mat'
model_CCSIM=reshape(Final,69,69,39,50);
%model_CCSIM=out1;
model_CCSIM(:,:,:,51)=out;

S=20*ones(1,51); S(51)=100;

dim = size(out);

abc=[1,2,3,4,5];

IX1=zeros(5,51);IX2=zeros(5,51);IX3=zeros(5,51);
for ii=1:5
    i=abc(ii);
    
    figure;
    ddd=squeeze(JS_dispat(i,:,:));
    [Y1_dispat,e_dispat] = cmdscale(double(ddd));
    x_MR=Y1_dispat(:,1); y_MR=Y1_dispat(:,2); z_MR=Y1_dispat(:,3);
    %%
    x_MR=x_MR-x_MR(51); y_MR=y_MR-y_MR(51); z_MR=z_MR-z_MR(51);
    %%
    [d1,IX1(ii,:)] = sort(ddd(51,:));
    
    %%
    scatter3(0,0,0,100,'k','filled');
    hold on;
    scatter3(x_MR,y_MR,z_MR, S, 'r','filled');
    
    hold on;
    ddd=squeeze(JS_CCSIM(i,:,:));
    [Y1_CCSIM,e_CCSIM] = cmdscale(double(ddd));
    x_CCSIM=Y1_CCSIM(:,1); y_CCSIM=Y1_CCSIM(:,2); z_CCSIM=Y1_CCSIM(:,3);
    %%
    x_CCSIM=x_CCSIM-x_CCSIM(51); y_CCSIM=y_CCSIM-y_CCSIM(51); z_CCSIM=z_CCSIM-z_CCSIM(51);
    [d2,IX2(ii,:)] = sort(ddd(51,:));
    
    scatter3(x_CCSIM,y_CCSIM,z_CCSIM, S,'b','filled');
    hold on;
%     ddd=squeeze(JS_sisim(i,:,:));
%     [Y1_sisim,e_sisim] = cmdscale(double(ddd));
%     x_sisim=Y1_sisim(:,1); y_sisim=Y1_sisim(:,2); z_sisim=Y1_sisim(:,3);
%     %%
%     x_sisim=x_sisim-x_sisim(51); y_sisim=y_sisim-y_sisim(51); z_sisim=z_sisim-z_sisim(51);
%     [d3,IX3(ii,:)] = sort(ddd(51,:));
%     
%     %S(IX3(2))=60; S(IX3(21))=60; S(IX3(31))=60; S(IX3(41))=60; S(IX3(51))=60;
%     scatter3(x_sisim,y_sisim,z_sisim,S,'r','filled');
%     hold on;
    scatter3(0,0,0,100,'k','filled');
    legend('Training image','DISPAT','CCSIM');
    title(['Multi Resolution= ' int2str(i)]);
    
    hold off;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure;
    S=20*ones(1,51); S(51)=100;
    S(IX1(ii,2))=30; S(IX1(ii,21))=30; S(IX1(ii,31))=30; S(IX1(ii,41))=30; S(IX1(ii,51))=30;
    
    
    
    scatter(x_MR,y_MR,S,'r','filled');
%     text(x_MR(IX1(ii,2)),y_MR(IX1(ii,2)),'\leftarrow 1','FontSize',10,'Color','g');
%     text(x_MR(IX1(ii,21)),y_MR(IX1(ii,21)),'\leftarrow 20','FontSize',10,'Color','g');
%     text(x_MR(IX1(ii,31)),y_MR(IX1(ii,31)),'\leftarrow 30','FontSize',10,'Color','g');
%     text(x_MR(IX1(ii,41)),y_MR(IX1(ii,41)),'\leftarrow 40','FontSize',10,'Color','g');
%     text(x_MR(IX1(ii,51)),y_MR(IX1(ii,51)),'\leftarrow 50','FontSize',10,'Color','g');
    hold on;
    S=20*ones(1,51); S(51)=100;
    S(IX2(ii,2))=30; S(IX2(ii,21))=30; S(IX2(ii,31))=30; S(IX2(ii,41))=30; S(IX2(ii,51))=30;
    scatter(x_CCSIM,y_CCSIM,S,'b','filled');
%     text(x_CCSIM(IX2(ii,2)),y_CCSIM(IX2(ii,2)),'\rightarrow 1','FontSize',10,'Color','b');
%     text(x_CCSIM(IX2(ii,21)),y_CCSIM(IX2(ii,21)),'\rightarrow 20','FontSize',10,'Color','b');
%     text(x_CCSIM(IX2(ii,31)),y_CCSIM(IX2(ii,31)),'\rightarrow 30','FontSize',10,'Color','b');
%     text(x_CCSIM(IX2(ii,41)),y_CCSIM(IX2(ii,41)),'\rightarrow 40','FontSize',10,'Color','b');
%     text(x_CCSIM(IX2(ii,51)),y_CCSIM(IX2(ii,51)),'\rightarrow 50','FontSize',10,'Color','b');
%     S=20*ones(1,51); S(51)=100;
%     S(IX3(ii,2))=30; S(IX3(ii,21))=30; S(IX3(ii,31))=30; S(IX3(ii,41))=30; S(IX3(ii,51))=30;
%     scatter(x_sisim,y_sisim,S,'r','filled');
%     text(x_sisim(IX3(ii,2)),y_sisim(IX3(ii,2)),'\leftarrow 1','FontSize',10,'Color','r');
%     text(x_sisim(IX3(ii,21)),y_sisim(IX3(ii,21)),'\leftarrow 20','FontSize',10,'Color','r');
%     text(x_sisim(IX3(ii,31)),y_sisim(IX3(ii,31)),'\leftarrow 30','FontSize',10,'Color','r');
%     text(x_sisim(IX3(ii,41)),y_sisim(IX3(ii,41)),'\leftarrow 40','FontSize',10,'Color','r');
%     text(x_sisim(IX3(ii,51)),y_sisim(IX3(ii,51)),'\leftarrow 50','FontSize',10,'Color','r');
    hold on;
    scatter(0,0,100,'k','filled');
    
    title(['x-y Plot & Multi Scale = ' int2str(i)]);
    legend('Dispat','CCSIM','training image');
    
    hold off;

%% find the corresponding models in the scatter plots


out_c=model_dispat(:,:,:,IX1(ii,1));
    out1 = my_resize(out_c,ceil(dim(1)/i),ceil(dim(2)/i),ceil(dim(3)/i)); 
    level = graythresh(out1);
    out1 = im2bw_3D(out1, level);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% use SGEMS to show these models
figure;
subplot(2,3,1)
%imshow(out1);
%imshow(out1(:,:,ceil(size(out1,3)/2)));
imshow(out1(:,:,ceil(size(out1,3))));
%h = vol3d('cdata',out1,'texture','2D');
%view(3); 
% % Update view since 'texture' = '2D'
%vol3d(h);  
%alphamap('rampdown'), alphamap('decrease'), alphamap('decrease')


title(['Ti, (0, 0), dispat, MR=' int2str(i)]);
subplot(2,3,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_c=model_dispat(:,:,:,IX1(ii,2));
    out1 = my_resize(out_c,ceil(dim(1)/i),ceil(dim(2)/i),ceil(dim(3)/i));  
    level = graythresh(out1);
    out1 = im2bw_3D(out1, level);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imshow(out1(:,:,ceil(size(out1,3)/2)));
imshow(out1(:,:,ceil(size(out1,3))));
%title(['closest to Ti, (' num2str(x_MR(IX1(2)),2) ',' num2str(y_MR(IX1(2)),2) ')'] );
title('closest to Ti');
for kk=1:4
    subplot(2,3,kk+2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_c=model_dispat(:,:,:,IX1(ii,11+10*kk));
    out1 = my_resize(out_c,ceil(dim(1)/i),ceil(dim(2)/i),ceil(dim(3)/i));   
    level = graythresh(out1);
    out1 = im2bw_3D(out1, level);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %imshow(out1(:,:,ceil(size(out1,3)/2)));
    imshow(out1(:,:,ceil(size(out1,3))));
    %title([int2str(11+10*kk-1) 'th closest to Ti, (' num2str(x_MR(IX1(11+10*kk)),2) ',' num2str(y_MR(IX1(11+10*kk)),2) ')'] );
    title([int2str(11+10*kk-1) 'th closest to Ti']);
end

%%
figure;
subplot(2,3,1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_c=model_CCSIM(:,:,:,IX2(ii,1));
    out1 = my_resize(out_c,ceil(dim(1)/i),ceil(dim(2)/i),ceil(dim(3)/i)); 
    level = graythresh(out1);
    out1 = im2bw_3D(out1, level);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imshow(out1(:,:,ceil(size(out1,3)/2)));
imshow(out1(:,:,ceil(size(out1,3))));
title(['Ti, (0, 0), CCSIM, MR=' int2str(i)]);
subplot(2,3,2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_c=model_CCSIM(:,:,:,IX2(ii,2));
    out1 = my_resize(out_c,ceil(dim(1)/i),ceil(dim(2)/i),ceil(dim(3)/i));
    level = graythresh(out1);
    out1 = im2bw_3D(out1, level);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%imshow(out1(:,:,ceil(size(out1,3)/2)));
imshow(out1(:,:,ceil(size(out1,3))));
%title(['closest to Ti, (' num2str(x_CCSIM(IX2(2)),2) ',' num2str(y_CCSIM(IX2(2)),2) ')'] );
title('closest to Ti');
for kk=1:4
    subplot(2,3,kk+2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
out_c=model_CCSIM(:,:,:,IX2(ii,11+10*kk));
    out1 = my_resize(out_c,ceil(dim(1)/i),ceil(dim(2)/i),ceil(dim(3)/i)); 
    level = graythresh(out1);
    out1 = im2bw_3D(out1, level);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %imshow(out1(:,:,ceil(size(out1,3)/2)));
    imshow(out1(:,:,ceil(size(out1,3))));
    %title([int2str(11+10*kk-1) 'th closest to Ti, (' num2str(x_CCSIM(IX2(11+10*kk)),2) ',' num2str(y_CCSIM(IX2(11+10*kk)),2) ')'] );
    title([int2str(11+10*kk-1) 'th closest to Ti']);
end
%%
% figure;
% subplot(2,3,1)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_c=model_sisim(:,:,:,IX3(ii,1));
%     out1 = my_resize(out_c,ceil(dim(1)/i),ceil(dim(2)/i),ceil(dim(3)/i));
%     level = graythresh(out1);
%     out1 = im2bw_3D(out1, level);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %imshow(out1(:,:,ceil(size(out1,3)/2)));
% imshow(out1(:,:,ceil(size(out1,3))));
% title(['Ti, (0, 0), snesim, MR=' int2str(i)]);
% subplot(2,3,2)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_c=model_sisim(:,:,:,IX3(ii,2));
%     out1 = my_resize(out_c,ceil(dim(1)/i),ceil(dim(2)/i),ceil(dim(3)/i)); 
%     level = graythresh(out1);
%     out1 = im2bw_3D(out1, level);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% imshow(out1(:,:,ceil(size(out1,3)/2)));
% %imshow(out1(:,:,ceil(size(out1,3))));
% %title(['closest to Ti, (' num2str(x_sisim(IX3(2)),2) ',' num2str(y_sisim(IX3(2)),2) ')'] );
% title('closest to Ti');
% for kk=1:4
%     subplot(2,3,kk+2)
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% out_c=model_sisim(:,:,:,IX3(ii,11+10*kk));
%     out1 = my_resize(out_c,ceil(dim(1)/i),ceil(dim(2)/i),ceil(dim(3)/i)); 
%     level = graythresh(out1);
%     out1 = im2bw_3D(out1, level);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %imshow(out1(:,:,ceil(size(out1,3)/2)));
%     imshow(out1(:,:,ceil(size(out1,3))));
%     %title([int2str(11+10*kk-1) 'th closest to Ti, (' num2str(x_sisim(IX3(11+10*kk)),2) ',' num2str(y_sisim(IX3(11+10*kk)),2) ')'] );
%     title([int2str(11+10*kk-1) 'th closest to Ti']);
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
S=20*ones(1,51); S(51)=100;
    S(IX1(ii,2))=30; S(IX1(ii,21))=30; S(IX1(ii,31))=30; S(IX1(ii,41))=30; S(IX1(ii,51))=30;
scatter(x_MR,z_MR,S,'r','filled');
% text(x_MR(IX1(ii,2)),z_MR(IX1(ii,2)),'\leftarrow 1','FontSize',10,'Color','g');
% text(x_MR(IX1(ii,21)),z_MR(IX1(ii,21)),'\leftarrow 20','FontSize',10,'Color','g');
% text(x_MR(IX1(ii,31)),z_MR(IX1(ii,31)),'\leftarrow 30','FontSize',10,'Color','g');
% text(x_MR(IX1(ii,41)),z_MR(IX1(ii,41)),'\leftarrow 40','FontSize',10,'Color','g');
% text(x_MR(IX1(ii,51)),z_MR(IX1(ii,51)),'\leftarrow 50','FontSize',10,'Color','g');
hold on;
S=20*ones(1,51); S(51)=100;
    S(IX2(ii,2))=30; S(IX2(ii,21))=30; S(IX2(ii,31))=30; S(IX2(ii,41))=30; S(IX2(ii,51))=30;
scatter(x_CCSIM,z_CCSIM,S,'b','filled');
% text(x_CCSIM(IX2(ii,2)),z_CCSIM(IX2(ii,2)),'\rightarrow 1','FontSize',10,'Color','b');
%     text(x_CCSIM(IX2(ii,21)),z_CCSIM(IX2(ii,21)),'\rightarrow 20','FontSize',10,'Color','b');
%     text(x_CCSIM(IX2(ii,31)),z_CCSIM(IX2(ii,31)),'\rightarrow 30','FontSize',10,'Color','b');
%     text(x_CCSIM(IX2(ii,41)),z_CCSIM(IX2(ii,41)),'\rightarrow 40','FontSize',10,'Color','b');
%     text(x_CCSIM(IX2(ii,51)),z_CCSIM(IX2(ii,51)),'\rightarrow 50','FontSize',10,'Color','b');
%  S=20*ones(1,51); S(51)=100;
%  S(IX3(ii,2))=30; S(IX3(ii,21))=30; S(IX3(ii,31))=30; S(IX3(ii,41))=30; S(IX3(ii,51))=30;
%  scatter(x_sisim,z_sisim,S,'r','filled');
%  text(x_sisim(IX3(ii,2)),z_sisim(IX3(ii,2)),'\leftarrow 1','FontSize',10,'Color','r');
%      text(x_sisim(IX3(ii,21)),z_sisim(IX3(ii,21)),'\leftarrow 20','FontSize',10,'Color','r');
%      text(x_sisim(IX3(ii,31)),z_sisim(IX3(ii,31)),'\leftarrow 30','FontSize',10,'Color','r');
%      text(x_sisim(IX3(ii,41)),z_sisim(IX3(ii,41)),'\leftarrow 40','FontSize',10,'Color','r');
%      text(x_sisim(IX3(ii,51)),z_sisim(IX3(ii,51)),'\leftarrow 50','FontSize',10,'Color','r');
title(['x-z Plot & Multi Scale = ' int2str(i)]);
hold on;
scatter(0,0,100,'k','filled');
legend('Dispat','CCSIM','training image');
hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
S=20*ones(1,51); S(51)=100;
S(IX1(ii,2))=30; S(IX1(ii,21))=30; S(IX1(ii,31))=30; S(IX1(ii,41))=30; S(IX1(ii,51))=30;
scatter(y_MR,z_MR,S,'r','filled');
% text(y_MR(IX1(ii,2)),z_MR(IX1(ii,2)),'\leftarrow 1','FontSize',10,'Color','g');
%     text(y_MR(IX1(ii,21)),z_MR(IX1(ii,21)),'\leftarrow 20','FontSize',10,'Color','g');
%     text(y_MR(IX1(ii,31)),z_MR(IX1(ii,31)),'\leftarrow 30','FontSize',10,'Color','g');
%     text(y_MR(IX1(ii,41)),z_MR(IX1(ii,41)),'\leftarrow 40','FontSize',10,'Color','g');
%     text(y_MR(IX1(ii,51)),z_MR(IX1(ii,51)),'\leftarrow 50','FontSize',10,'Color','g');
hold on;
S=20*ones(1,51); S(51)=100;
    S(IX2(ii,2))=30; S(IX2(ii,21))=30; S(IX2(ii,31))=30; S(IX2(ii,41))=30; S(IX2(ii,51))=30;
scatter(y_CCSIM,z_CCSIM,S,'b','filled');
% text(y_CCSIM(IX2(ii,2)),z_CCSIM(IX2(ii,2)),'\rightarrow 1','FontSize',10,'Color','b');
%     text(y_CCSIM(IX2(ii,21)),z_CCSIM(IX2(ii,21)),'\rightarrow 20','FontSize',10,'Color','b');
%     text(y_CCSIM(IX2(ii,31)),z_CCSIM(IX2(ii,31)),'\rightarrow 30','FontSize',10,'Color','b');
%     text(y_CCSIM(IX2(ii,41)),z_CCSIM(IX2(ii,41)),'\rightarrow 40','FontSize',10,'Color','b');
%     text(y_CCSIM(IX2(ii,51)),z_CCSIM(IX2(ii,51)),'\rightarrow 50','FontSize',10,'Color','b');
%  S=20*ones(1,51); S(51)=100;
%      S(IX3(ii,2))=30; S(IX3(ii,21))=30; S(IX3(ii,31))=30; S(IX3(ii,41))=30; S(IX3(ii,51))=30;
%  scatter(y_sisim,z_sisim,S,'g','filled');
%  text(y_sisim(IX3(ii,2)),z_sisim(IX3(ii,2)),'\rightarrow 1','FontSize',10,'Color','r');
%      text(y_sisim(IX3(ii,21)),z_sisim(IX3(ii,21)),'\rightarrow 20','FontSize',10,'Color','r');
%      text(y_sisim(IX3(ii,31)),z_sisim(IX3(ii,31)),'\rightarrow 30','FontSize',10,'Color','r');
%      text(y_sisim(IX3(ii,41)),z_sisim(IX3(ii,41)),'\rightarrow 40','FontSize',10,'Color','r');
%      text(y_sisim(IX3(ii,51)),z_sisim(IX3(ii,51)),'\rightarrow 50','FontSize',10,'Color','r');

title(['y-z Plot & Multi Scale = ' int2str(i)]);
hold on;
scatter(0,0,100,'k','filled');
legend('Dispat','CCSIM','training image');
hold off; 

size_dispatE=size(Y1_dispat,2);
size_CCSIME=size(Y1_CCSIM,2);
% size_sisimE=size(Y1_sisim,2);
for k=1:size_dispatE
    e_dispat_new(k)=sum(e_dispat(1:k))/sum(e_dispat(1:size_dispatE)); 
end
for k=1:size_CCSIME
    e_CCSIM_new(k)=sum(e_CCSIM(1:k))/sum(e_CCSIM(1:size_CCSIME));
end
%  for k=1:size_sisimE
%      e_sisim_new(k)=sum(e_sisim(1:k))/sum(e_sisim(1:size_sisimE));
%  end


figure;
plot(e_dispat_new,'k');
%text(3,e_dispat_new(3),num2str(e_dispat_new(3)));
hold on;
plot(e_CCSIM_new,'b');
%text(3,e_CCSIM_new(3),num2str(e_CCSIM_new(3)));
%  plot(e_sisim_new,'r');
%text(3,e_sisim_new(3),num2str(e_sisim_new(3)));
xlabel('Dimension');
ylabel('Cumulative');
legend(['Dispat ' num2str(e_dispat_new(3))],['CCSIM ' num2str(e_CCSIM_new(3))],'Location', 'SouthEast');
title(['Multi Resolution= ' int2str(i)]);
hold off;

end



%save('ranking_KL.mat', 'IX1','IX2','IX3');
