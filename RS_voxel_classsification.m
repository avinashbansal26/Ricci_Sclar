%% For  0 to 5 fiber image
%  [All_Tensor_Coff]=Fibers_specific_iteration();
 [All_Tensor_Coff]=Fibers_0_to_5_specific_iteration();
%    [All_Tensor_Coff]=FibersSimulation();
%  [All_Tensor_Coff]=geodesic_simulation_cross();
%   [All_Tensor_Coff]=cross_fiber();
UnitVectors ;
n_grad=81;
GradientOrientations=g(1:n_grad,:);
flag=1;
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
        FinSum=zeros(3,3);
        for ii=1:size(v,1)
            if (v(ii,1)>=0)  % Choosing all direction on hemisphere
                count=count+1;
                p_dirn(count,:)=v(ii,:);  % only postive dirn / dirn on hemisphere only
                F(:,:,count)=double(F_metric (TD(1),TD(2),TD(3),TD(4),TD(5),TD(6),TD(7),TD(8),TD(9),TD(10),TD(11),...
                    TD(12),TD(13), TD(14), TD(15), v(ii,1),v(ii,2),v(ii,3)));
                Prin_Max_Eig(count)= max(eig((squeeze(F(:,:,count)))));
                Fin(:,:,i,j,count)=squeeze(F(:,:,count));
                FinSum=FinSum+squeeze(F(:,:,count));
            end
        end
        FSumAll(:,:,i,j)=FinSum;
    end
end
% figure; plotDTI(FSumAll,.5);

% for ii=1:max_num
%     figure; plotDTI(squeeze(Fin(:,:,:,:,ii)),2);
%% Mertic tensor g and gi (g inverse)
% g=squeeze(Fin(:,:,:,:,ii));

figure;
for ii=1:8:size(p_dirn,1)
    g=squeeze(Fin(:,:,:,:,ii));
    plotDTI(g,3);
end
title('All Finsler Metrix');

Fsum=zeros(3,3,size(D,2),size(D,3));
count=0;
for ii=1:8:size(p_dirn,1)
%     if flag==1
%         g=FSumAll;
%         flag=0;
%     else
count=count+1;
        g=squeeze(Fin(:,:,:,:,ii));
        Fsum=Fsum+g;
        figure; plotDTI(g,3);
        title(ii)
        for i=1:size(g,3)
            for j=1:size(g,4)
                gi(:,:,i,j)=inv(g(:,:,i,j));
            end
        end
        
        
        
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
        title(ii);
        % figure; imagesc(RS);
        RS_All(1:size(gi_inter,3),1:size(gi_inter,4),count)=RS_by_RicTensor;
%     end
end
RS_Sum=zeros(size(gi_inter,3),size(gi_inter,4));
figure; plotDTI(Fsum,20);
title('FSum DTI');
for ii=1:size(RS_All,3)
    RS_Sum=RS_Sum + squeeze(RS_All(:,:,ii));
end
figure; surf(RS_Sum);
title('RS Sum Surface');

cluster=3;
idx = kmeans(RS_Sum,cluster);
figure; gscatter(RS_Sum,idx);

%% K means Algorithm
count=0;
for i=1:size(RS_All,1)
    for j=1:size(RS_All,2)
        count=count+1;
        Kmean_RS(count,:)= squeeze(RS_All(i,j,:));
    end
end
cluster=3;
idx = kmeans(Kmean_RS,cluster);
figure; gscatter(Kmean_RS,idx);
title('Sample Data');

count=0;
for i=1:size(RS_All,1)
    for j=1:size(RS_All,2)
        count=count+1;
        RS_Data(i,j)= idx(count);
    end
end
figure;imagesc(RS_Data);

