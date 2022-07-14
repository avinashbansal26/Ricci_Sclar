function [TC_region]=open_HARDICeitec()

% [V,info]=ReadData3D('C:\Users\HP\Desktop\dwmriDATA\Diffusion_dir64_b1500\Diffusion\2621B_03d_ep2d_diff_mddw_2iso_dir64_b1500_AP.nii') ;%'340B_data.nii' 4Ddwi_b1000  ceitec.nii
%  [V,info]=ReadData3D('C:\Users\kaushik\AllCodes\Codes_Fourth_order\cc.nii'); %'340B_data.nii' 4Ddwi_b1000  ceitec.nii
[V]=ReadData3D('C:\Users\bansala\Desktop\Avinash\Diffusion\2621B_03d_ep2d_diff_mddw_2iso_dir64_b1500_AP.nii');

% V : The 3D Volume
% info : Struct with info about the data
% Always the following fields are present
% info.Filename : Name of file
% info.Dimensions : Dimensions of Volume
% info.PixelDimensions : Size of one pixel / voxel

UnitVectorsHardiCeitec;
% phantom
n_grad=74;%64 %73
f=sqrt(2);
lsize=16;
GradientOrientations= g1([1:n_grad],:);
%  pause;
bvalue=1000;
%   b_value=[0;ones(n_grad-9,1)*1500;0;0;0;0;0;0;0;0;0];
b_value= [0;ones(n_grad-10,1)*1500;0;0;0;0;0;0;0;0;0];
%  S=ones(size(V,1),size(V,2),size(V,3),size(GradientOrientations,1));
order=4;

S=V;
size(V);
order=4;
x_ly_no=51;
y_ly_no=56;
z_ly_no=36;
xsize=1:114;
ysize=1:114;
zsize=1:70;
tS=squeeze(S(:,y_ly_no,:,1));
pause;
%   tS=imrotate(tS,90);
figure;imagesc(imrotate(tS,90));
figure;imagesc(imrotate(tS(xsize,zsize),90));
% size(S)
%   pause;
G=constructMatrixOfMonomials(GradientOrientations, order); %computes G from section 5.1 (ISBI'10)
%Construct set of polynomial coefficients C
C=constructSetOf321Polynomials(order)'; %computes C from section 5.1 (ISBI'10)
P=G*C;
% pause;
P=[-diag(b_value)*P ones(size(GradientOrientations,1),1)];
BG=constructMatrixOfIntegrals(GradientOrientations, order, 100);
B=BG*C;


for i=1:size(S,1) %x_ly_no:x_ly_no
    for j=y_ly_no:y_ly_no
        for k=1:size(S,3)
            for gcounter=1:size(S,4)
                if(S(i,j,k,gcounter)==0)
                    S(i,j,k,gcounter)=0.0001;
                end
            end
            y1=(squeeze(log(double(S(i,j,k,:)))));
            if(y1=='Inf' )
                y1(:)=9999;
            end
            if(y1==0 )
                y1(:)=0.0001;
            end
            x1=lsqnonneg(P, y1);
            if(x1=='Inf' )
                x1(:)=9999;
            end
            if(x1==0 )
                x1(:)=0.0001;
            end
            CT_Coeff = C * x1([1:321]);
            TC(:,i,k)=CT_Coeff;
        end
    end
end
disp('hi');
figure; plotTensors(squeeze(TC(:,xsize,zsize)),1,[321  1]);
TC_region=TC(:,xsize,zsize);
  TC_region=TC(:,40:80,20:50);
%  figure; plotTensors(TC_region,1,[321  1]);
end
