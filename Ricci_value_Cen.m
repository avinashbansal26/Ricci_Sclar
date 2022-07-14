%% For  0 to 5 fiber image

[All_Tensor_Coff]=FibersSimulation();
for i=1:size(All_Tensor_Coff,2)+4
    for j=1:size(All_Tensor_Coff,3)+4
      Tensor_Coff(:,i,j)=  1.0e-05 * [0.1016   -0.0000    0.2031    0.0000    0.1016   -0.0000 ...
          -0.0000    0.0000   -0.0000    0.2030    0.0000  0.2030   -0.0000    0.0000    0.1016];
    end
end
Tensor_Coff(:,3:size(All_Tensor_Coff,2)+2,3:size(All_Tensor_Coff,3)+2)=All_Tensor_Coff;

UnitVectors ;
n_grad=81;
GradientOrientations=[1 0 0; g(1:n_grad,:)];

%% Next set of lines for synthetic image (Another simmulation crossind in 3 fibers Circle based)
%    [All_Tensor_Coff]=geodesic_simulation_cross();
%    UnitVectors;
%    n_grad=81;
%    GradientOrientations=[1 0 0;g([1:n_grad],:)];

v= GradientOrientations;  % 64 gradient direction
syms x11 x12 x13 x22 x23 x33;
syms xL11 xL12 xL13 xL22 xL23 xL33;
syms xR11 xR12 xR13 xR22 xR23 xR33;
syms xU11 xU12 xU13 xU22 xU23 xU33;
syms xD11 xD12 xD13 xD22 xD23 xD33;
syms xLL11 xLL12 xLL13 xLL22 xLL23 xLL33;
syms xLU11 xLU12 xLU13 xLU22 xLU23 xLU33;
syms xLD11 xLD12 xLD13 xLD22 xLD23 xLD33;
syms xRR11 xRR12 xRR13 xRR22 xRR23 xRR33;
syms xRU11 xRU12 xRU13 xRU22 xRU23 xRU33;
syms xRD11 xRD12 xRD13 xRD22 xRD23 xRD33;
syms xUU11 xUU12 xUU13 xUU22 xUU23 xUU33;
syms xDD11 xDD12 xDD13 xDD22 xDD23 xDD33;

%  load('RS_Cen.mat');
RS_Cen(x11,x12,x13,x22,x23,x33,xL11,xL12,xL13,xL22,xL23,xL33,xR11,xR12,xR13,xR22,xR23,xR33,...
 xU11,xU12,xU13,xU22,xU23,xU33,xD11,xD12,xD13,xD22,xD23,xD33,xLL11,xLL12,xLL13,xLL22,xLL23,xLL33,...
 xLU11,xLU12,xLU13,xLU22,xLU23,xLU33,xLD11,xLD12,xLD13,xLD22,xLD23,xLD33,xRR11,xRR12,xRR13,xRR22,xRR23,xRR33,...
 xRU11,xRU12,xRU13,xRU22,xRU23,xRU33,xRD11,xRD12,xRD13,xRD22,xRD23,xRD33,xUU11,xUU12,xUU13,xUU22,xUU23,xUU33,...
 xDD11,xDD12,xDD13,xDD22,xDD23,xDD33)= Ricci_Scalar_NM_Cen_var();
pause;
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
        max_num=3;
        [maxi,index]=maxk(Prin_Max_Eig,max_num);
%         FSum=zeros(3,3);
%         Fin=zeros(3,3,max_num);
        for ii=1:max_num
            Fin(:,:,i,j,ii)=(squeeze(F(:,:,index(ii))));
%             FSum=FSum+squeeze(F(:,:,ii));
        end
%         FSumAll(:,:,i,j)=FSum/count;         % Mean of Finsler/ diffusion matrix
    end
end
figure;
for ii=1:max_num   
      plotDTI(squeeze(Fin(:,:,3:size(D,2)-2,3:size(D,3)-2,ii)),12);    
end

 FSumAll=zeros(3,3,size(D,2),size(D,3));
 for ii=1:max_num
     FSumAll= squeeze(Fin(:,:,1:size(D,2),1:size(D,3),ii));
%   FSumAll=  FSumAll+squeeze(Fin(:,:,1:size(D,2),1:size(D,3),ii));
% end 

% FA_F=GFA(FSumAll);
% figure; plotDTI(FSumAll,20); 
for x=3:size(D,2)-2
    for y=3:size(D,3)-2
        i=x; j=y;
        [x11, x12, x13, x22, x23, x33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x+1; j=y;
        [xR11, xR12, xR13, xR22, xR23, xR33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x-1; j=y;
        [xL11, xL12, xL13, xL22, xL23, xL33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x; j=y+1;
        [xU11, xU12, xU13, xU22, xU23, xU33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x; j=y-1;
        [xD11, xD12, xD13, xD22, xD23, xD33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x+1; j=y+1;
        [xRU11, xRU12, xRU13, xRU22, xRU23, xRU33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x+1; j=y-1;
        [xRD11, xRD12, xRD13, xRD22, xRD23, xRD33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x-1; j=y+1;
        [xLU11, xLU12, xLU13, xLU22, xLU23, xLU33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
         i=x-1; j=y-1;
        [xLD11, xLD12, xLD13, xLD22, xLD23, xLD33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x+2; j=y;
        [xRR11, xRR12, xRR13, xRR22, xRR23, xRR33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x-2; j=y;
        [xLL11, xLL12, xLL13, xLL22, xLL23, xLL33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        i=x; j=y+2;
        [xUU11, xUU12, xUU13, xUU22, xUU23, xUU33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));  
        i=x; j=y-2;
        [xDD11, xDD12, xDD13, xDD22, xDD23, xDD33]=deal(FSumAll(1,1,i,j),FSumAll(1,2,i,j),FSumAll(1,3,i,j),...
            FSumAll(2,2,i,j),FSumAll(2,3,i,j),FSumAll(3,3,i,j));
        
        Ricci_Sclar(x,y)=double(RS_Cen(x11,x12,x13,x22,x23,x33,xL11,xL12,xL13,xL22,xL23,xL33,xR11,xR12,xR13,xR22,xR23,xR33,...
                    xU11,xU12,xU13,xU22,xU23,xU33,xD11,xD12,xD13,xD22,xD23,xD33,xLL11,xLL12,xLL13,xLL22,xLL23,xLL33,...
                    xLU11,xLU12,xLU13,xLU22,xLU23,xLU33,xLD11,xLD12,xLD13,xLD22,xLD23,xLD33,xRR11,xRR12,xRR13,xRR22,xRR23,xRR33,...
                    xRU11,xRU12,xRU13,xRU22,xRU23,xRU33,xRD11,xRD12,xRD13,xRD22,xRD23,xRD33,xUU11,xUU12,xUU13,xUU22,xUU23,xUU33,...
                    xDD11,xDD12,xDD13,xDD22,xDD23,xDD33));
        if Ricci_Sclar(x,y)<0
            rgbImage(x, y, :) = [1, 0, 0]; %Red
        elseif Ricci_Sclar(x,y)>0
            rgbImage(x, y, :) = [0, 1, 0]; %Green
        else
            rgbImage(x, y, :) = [0, 0, 1]; %Blue
        end           
    end
end
Ric_Sclr(:,:,ii)=Ricci_Sclar(3:size(D,2)-2, 3:size(D,3)-2);
 end
  rgbImage = rgbImage(3:size(D,2)-2, 3:size(D,3)-2, :);
Ric_Sclr_Sum=zeros(size(D,2)-4, size(D,3)-4);
for ii=1:max_num
    Ric_Sclr_Sum =Ric_Sclr_Sum + squeeze(Ric_Sclr(:,:,ii));
figure; surf(squeeze(Ric_Sclr(:,:,ii)));
end
% figure; surf(Ricci_Sclar);
%  figure; imagesc(Ricci_Sclar);
% figure; imagesc(rgbImage);
figure; surf(Ric_Sclr_Sum);
 figure; imagesc(Ric_Sclr_Sum);
 
Ric_Sclr_Sum
figure; imagesc(rgbImage);


