%% For  0 to 5 fiber image
[All_Tensor_Coff]=Fibers_specific_iteration();
%    [All_Tensor_Coff]=Fibers_0_to_5_specific_iteration();
%      [All_Tensor_Coff]=FibersSimulation();
%    [All_Tensor_Coff]=geodesic_simulation_cross();
%       [All_Tensor_Coff]=cross_fiber();
UnitVectors ;
n_grad=81;
GradientOrientations=[g(1:n_grad,:)];
figure; plot3(GradientOrientations(:,1),GradientOrientations(:,2),GradientOrientations(:,3),'o');
title('Direction for Syntheatic Images');

%% For real brain data uncomment this section
%{
[All_Tensor_Coff]= open_HARDICeitec();
UnitVectorsHardiCeitec;
n_grad=65;%64 %73
GradientOrientations= g1([2:n_grad],:);
figure; plot3(GradientOrientations(:,1),GradientOrientations(:,2),GradientOrientations(:,3),'o');
title('Direction for Real Images');
%}

%% Next set of lines for synthetic image (Another simmulation crossind in 3 fibers Circle based)
%    [All_Tensor_Coff]=geodesic_simulation_cross();
%    UnitVectors;
%    n_grad=81;
%    GradientOrientations=[1 0 0;g([1:n_grad],:)];

v= GradientOrientations;  % 64 gradient direction
load('F_metric.mat');
% F_metric=Finsler_metric_as_variables();  % symbolic matrix which needs 15 cofficent and [1 3] vector

D=All_Tensor_Coff; % It will run for full 2D brain image and time consuling.
for i=1:size(D,2)
    for j=1:size(D,3)
        %% code for Finsler Fractional Anisotropy (FFA)
        TD=D(:,i,j)+.0001; % Adding small value(0.0001) to avoid division by zero error
        count=0;
        for ii=1:size(v,1)
            if (v(ii,1)>=0)  % Choosing all direction on hemisphere
                count=count+1;
                F(:,:,count)=double(F_metric (TD(1),TD(2),TD(3),TD(4),TD(5),TD(6),TD(7),TD(8),TD(9),TD(10),TD(11),...
                    TD(12),TD(13), TD(14), TD(15), v(ii,1),v(ii,2),v(ii,3)));
                Prin_Max_Eig(count)= max(eig((squeeze(F(:,:,count)))));
            end
        end
        max_num=6;
        [maxi,index]=maxk(Prin_Max_Eig,max_num);
        FSum=zeros(3,3);
        for ii=1:max_num
            FSum=FSum+squeeze(F(:,:,index(ii)));
            Fin(:,:,i,j,ii)=squeeze(F(:,:,index(ii)));
        end
        FSumAll(:,:,i,j)=FSum/count;         % Mean of Finsler/ diffusion matrix
    end
end
figure; plotDTI(FSumAll,.5);
title('Finsler Sum');
figure;
for ii=1:max_num
    plotDTI(squeeze(Fin(:,:,:,:,ii)),12);
end
title('All Fin overlap');
%% Mertic tensor g and gi (g inverse)
% g=squeeze(Fin(:,:,:,:,ii));
g=FSumAll;
%% DTI data just for experiment.
%{
load('gDTI_Data.mat')
g=squeeze(gDTI(:,:,40:90,40:100,29)+diag([.001,.001,.001]));
figure; plotDTI(g,.002);
title('only DTI Data ');
FA_gDTI=GFA(g);
figure; imagesc(imrotate(FA_gDTI,90));
title('Fractinal Anisotropy of DTI Data ');
%}
for i=1:size(g,3)
    for j=1:size(g,4)
        gi(:,:,i,j)=inv(g(:,:,i,j));
        giwis(i,j)=sum(sum( squeeze(gi(:,:,i,j))));
    end
end
figure; imagesc(imrotate(giwis,90));
title('g inverse without interpolation');
FA_FsumAll=GFA(FSumAll);
figure; imagesc(imrotate(FA_FsumAll,90));
title('FA FsumAll');
figure; surf(FA_FsumAll);
title('FA FsumAll');


%% Interpolation of g, gi (g inverse)  and dgx (derivative of g along x (x=x1,x2,x3))
[X,Y]=meshgrid(1:size(g,4),1:size(g,3));
[Xq,Yq] = meshgrid(1:0.25:size(g,4), 1:0.25:size(g,3));
for ei=1:3
    for ej=1:3
        g_inter(ei,ej,:,:)=interp2(X,Y,squeeze(g(ei,ej,:,:)),Xq,Yq);
        gi_inter(ei,ej,:,:)=interp2(X,Y,squeeze(gi(ei,ej,:,:)),Xq,Yq);
    end
end

%% Central derivative of g along X and Y
for x=2:size(gi_inter,3)-1
    for y=2:size(gi_inter,4)-1
        for ei=1:3
            for ej=1:3
                dgx_inter(ei,ej,1,x,y)=(1/2)*(squeeze(g_inter(ei,ej,x+1,y))-squeeze(g_inter(ei,ej,x-1,y)));
                dgx_inter(ei,ej,2,x,y)=(1/2)*(squeeze(g_inter(ei,ej,x,y+1))-squeeze(g_inter(ei,ej,x,y-1)));
                dgx_inter(ei,ej,3,x,y)=0;
            end
        end
    end
end
cdgx_inter=zeros(3,3,3,size(gi_inter,3),size(gi_inter,4));
cdgx_inter(:,:,:,1:size(gi_inter,3)-1,1:size(gi_inter,4)-1)=dgx_inter;
lopsz=3;
%% 1st derivative_inter sclar
for x=1:size(gi_inter,3)
    for y=1:size(gi_inter,4)
        c=0; g=0; gi=0;
        for i=1:lopsz
            for j=1:lopsz
                g=g+g_inter(i,j,x,y);
                gi=gi+gi_inter(i,j,x,y);
                for k=1:lopsz
                    c=c+cdgx_inter(i,j,k,x,y);
                end
            end
        end
        ds(x,y)=c;
        gs(x,y)=g;
        gis(x,y)=gi;
    end
end
figure; imagesc(imrotate(gs,90));
title('g inter Symbol');
figure; surf(gs);
title('g inter Symbol');
figure; imagesc(imrotate(gis,90));
title('gi g inverse inter Symbol');
figure; surf(gis);
title('gi g inverse inter Symbol');
figure; imagesc(imrotate(ds,90));
title('1st derivative Symbol');
figure; surf(ds);
title('1st derivative Symbol');




%% Christoffel Symbol
lopsz=3;
for x=1:size(gi_inter,3)
    for y=1:size(gi_inter,4)
        for i=1:lopsz
            for j=1:lopsz
                for k=1:lopsz
                    c=0;
                    for l=1:lopsz
                        c = c + (1/2)*gi_inter(k,l,x,y)*(cdgx_inter(i,l,j,x,y)+cdgx_inter(j,l,i,x,y)-cdgx_inter(i,j,l,x,y));
                    end
                    crst_inter(i,j,k,x,y)=c;
                end
            end
        end
    end
end

%% crst_inter sclar

for x=1:size(gi_inter,3)
    for y=1:size(gi_inter,4)
        c=0;
        for i=1:lopsz
            for j=1:lopsz
                for k=1:lopsz
                    c=c+crst_inter(i,j,k,x,y);
                end
            end
        end
        cs(x,y)=c;
    end
end
figure; imagesc(imrotate(cs,90));
title('Christoffel Symbol');
figure; surf(cs);
title('Christoffel Symbol');

%% Central derivative of Christoffel Symbol along X and Y
for x=2:size(gi_inter,3)-1
    for y=2:size(gi_inter,4)-1
        for ei=1:3
            for ej=1:3
                for ek=1:3
                    temp_dcrstx_inter(ei,ej,ek,1,x,y)=(1/2)*(squeeze(crst_inter(ei,ej,ek,x+1,y))-squeeze(crst_inter(ei,ej,ek,x-1,y)));
                    temp_dcrstx_inter(ei,ej,ek,2,x,y)=(1/2)*(squeeze(crst_inter(ei,ej,ek,x,y+1))-squeeze(crst_inter(ei,ej,ek,x,y-1)));
                    temp_dcrstx_inter(ei,ej,ek,3,x,y)=0;
                end
            end
        end
    end
end
dcrstx_inter=zeros(3,3,3,3,size(gi_inter,3),size(gi_inter,4));
dcrstx_inter(:,:,:,:,1:size(gi_inter,3)-1,1:size(gi_inter,4)-1)=temp_dcrstx_inter;

%% Ricci Tensor
for x=1:size(gi_inter,3)
    for y=1:size(gi_inter,4)
        for i=1:lopsz
            for j=1:lopsz
                A=0;
                B=0;
                C=0;
                D=0;
                for a=1:lopsz
                    C=C+dcrstx_inter(i,j,a,a,x,y);
                    D=D+dcrstx_inter(a,j,a,i,x,y);
                    for b=1:lopsz
                        A=A+crst_inter(a,b,a,x,y)*crst_inter(i,j,b,x,y);
                        B=B+crst_inter(i,b,a,x,y)*crst_inter(a,j,b,x,y);
                    end
                end
                RicTensor_inter(i,j,x,y)=A-B+C-D;
            end
        end
    end
end
%% Ricci Scalar by Ricci Tensor
RS_by_RicTensor= (sum(sum(gi_inter.*RicTensor_inter))); % compute Ricci Scalar from Ricci Tensor
RS(1:size(gi_inter,3),1:size(gi_inter,4))=RS_by_RicTensor;
figure; surf(RS(3:size(gi_inter,3)-3,3:size(gi_inter,4)-3));
title('Surf of Ricci Sclar')
%  figure; imagesc(imrotate(RS,90));
% RS_All(1:size(gi_inter,3),1:size(gi_inter,4),ii)=RS_by_RicTensor;
% end
% RS_Sum=zeros(size(gi_inter,3),size(gi_inter,4));
% for ii=1:max_num
%     RS_Sum=RS_Sum+ squeeze(RS_All(:,:,ii));
% end
% figure; surf(RS_Sum);

%% RGB Image
for i=3:size(gi_inter,3)-3
    for j=3:size(gi_inter,4)-3
        if(RS(i,j)>0)
            RGB(i,j,:)=[1 0 0];
        elseif RS(i,j)==0
            RGB(i,j,:)=[0 1 0];
        elseif RS(i,j)<0
            RGB(i,j,:)=[0 0 1];
        end
    end
end

figure; imagesc(imrotate(RGB,90));
title('RGB image of Ricci Sclar')
