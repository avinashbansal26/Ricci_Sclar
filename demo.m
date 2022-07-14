function [g_inter]=demo()
T_X=2.5*diag([0.9511    0.2629    0.1625]);
T_Y=2.5*diag([0.2629  0.9511  0.1625]);
thx=45; thy=45;
rot_x=[1   0  0; 0   cos(thx)  -sin(thx); 0  sin(thx)  cos(thx)];
rot_y=[cos(thy)  0  sin(thy); 0 1 0; -sin(thy) 0 cos(thy)];
R=rot_y*rot_x;
T_45= R'*T_X*R;
thx=90; thy=90;
rot_x=[1   0  0; 0   cos(thx)  -sin(thx); 0  sin(thx)  cos(thx)];
rot_y=[cos(thy)  0  sin(thy); 0 1 0; -sin(thy) 0 cos(thy)];
R=rot_y*rot_x;
T_135=R'*T_45*R;
T_S= diag([ 1/sqrt(3) 1/sqrt(3)   1/sqrt(3)]);
T_XY=.5*(T_X+T_Y);
T_30=[4.9566 5.1116 0; 5.1116 8.1247 0;0 0 1.0000]/5;
T_45=[ 2 1 0; 1 1 0; 0 0 .5];
T_45_C=[ 2 15 0; 15 1 0; 0 0 .5]/15;
T_135=[ 2 -1 0; -1 1 0; 0 0 .5];
T_135_C=[ 2 -15 0; -15 1 0; 0 0 .5]/15;
T_150=[4.9566 -5.1116 0; -5.1116 8.1247 0;0 0 1.0000]/5;

lopsz=5;
%% demo 1
S=zeros(3,3,lopsz,lopsz);
for i=1:lopsz
    for j=1:lopsz    
        S(:,:,i,j)=T_S;
    end
end
figure; plotDTI(S);

%% demo 2
S=zeros(3,3,lopsz,lopsz);
for i=1:lopsz
    for j=1:lopsz    
        S(:,:,i,j)=T_45;
    end
end
figure; plotDTI(S,2);

%% demo 3
S=zeros(3,3,lopsz,lopsz);
for i=1:lopsz
    for j=1:lopsz    
        if (rem(i,2)==0 ||  rem(j,3)==0)
        S(:,:,i,j)=T_S;
        end
        if (rem(i,3)==0 ||  rem(j,5)==0)
        S(:,:,i,j)=2*T_S;
        end
        if (rem(i,4)==0 ||  rem(j,2)==0)
        S(:,:,i,j)=T_S/2;
        end
    end
end
figure; plotDTI(S,2);

% %% demo 4
% S=zeros(3,3,lopsz,lopsz);
% for i=1:lopsz
%     for j=1:lopsz    
%         if (rem(i,2)==0 ||  rem(j,3)==0)
%         S(:,:,i,j)=T_S;
%         end
%         if (rem(i,3)==0 ||  rem(j,5)==0)
%         S(:,:,i,j)=T_45/3;
%         end
%         if (rem(i,4)==0 ||  rem(j,2)==0)
%         S(:,:,i,j)=T_135/2;
%         end
%         if (rem(i,5)==0 ||  rem(j,4)==0)
%         S(:,:,i,j)=2*T_S;
%         end
%     end
% end
% figure; plotDTI(S,2);
% 
% %% demo 5
% S=zeros(3,3,lopsz,lopsz);
% S(:,:,1,1)=T_45; S(:,:,1,2)=T_30; S(:,:,1,3)=T_Y; S(:,:,1,4)=T_150; S(:,:,1,5)=T_135;
% S(:,:,2,1)=T_30; S(:,:,2,2)=T_XY; S(:,:,2,3)=T_XY; S(:,:,2,4)=T_XY; S(:,:,2,5)=T_150;
% S(:,:,3,1)=T_X;   S(:,:,3,2)=T_XY; S(:,:,3,3)=T_XY; S(:,:,3,4)=T_XY; S(:,:,3,5)=T_X;
% S(:,:,4,1)=T_150; S(:,:,4,2)=T_XY; S(:,:,4,3)=T_XY; S(:,:,4,4)=T_XY; S(:,:,4,5)=T_30;
% S(:,:,5,1)=T_135; S(:,:,5,2)=T_150; S(:,:,5,3)=T_Y; S(:,:,5,4)=T_30; S(:,:,5,5)=T_45;
% figure; plotDTI(S,2);
% 
% %% demo 6
% S=zeros(3,3,lopsz,lopsz);
% S(:,:,1,1)=T_Y; S(:,:,1,2)=T_150; S(:,:,1,3)=T_X; S(:,:,1,4)=T_30; S(:,:,1,5)=T_Y;
% S(:,:,2,1)=T_30; S(:,:,2,2)=T_XY; S(:,:,2,3)=T_XY; S(:,:,2,4)=T_XY; S(:,:,2,5)=T_150;
% S(:,:,3,1)=T_X;   S(:,:,3,2)=T_XY; S(:,:,3,3)=2*T_S; S(:,:,3,4)=T_XY; S(:,:,3,5)=T_X;
% S(:,:,4,1)=T_150; S(:,:,4,2)=T_XY; S(:,:,4,3)=T_XY; S(:,:,4,4)=T_XY; S(:,:,4,5)=T_30;
% S(:,:,5,1)=T_Y; S(:,:,5,2)=T_30; S(:,:,5,3)=T_X; S(:,:,5,4)=T_150; S(:,:,5,5)=T_Y;
% figure; plotDTI(S,2);



%% interpolation figure; plotDTI(g_inter,2)
g=S;
inter_size=0.25;
[X,Y]=meshgrid(1:size(g,4),1:size(g,3));
[Xq,Yq] = meshgrid(1:inter_size:size(g,4), 1:inter_size:size(g,3));
for ei=1:3
    for ej=1:3
        g_inter(ei,ej,:,:)=interp2(X,Y,squeeze(g(ei,ej,:,:)),Xq,Yq); 
    end
end
 figure; plotDTI(g_inter,2)
end