
% for ii=1:max_num
%     figure; plotDTI(squeeze(Fin(:,:,:,:,ii)),2);
%% Mertic tensor g and gi (g inverse)
% g=squeeze(Fin(:,:,:,:,ii));
g=demo;
for i=1:size(g,3)
    for j=1:size(g,4)
        gi_inter(:,:,i,j)=inv(g(:,:,i,j));
    end
end

%% Central derivative of g along X and Y
for x=2:size(g,3)-1
    for y=2:size(g,4)-1
        for ei=1:3
            for ej=1:3
                dgx_inter(ei,ej,1,x,y)=(1/2)*(squeeze(g(ei,ej,x+1,y))-squeeze(g(ei,ej,x-1,y)));
                dgx_inter(ei,ej,2,x,y)=(1/2)*(squeeze(g(ei,ej,x,y+1))-squeeze(g(ei,ej,x,y-1)));
                dgx_inter(ei,ej,3,x,y)=0;
            end
        end
    end
end
% cdgx=zeros(3,3,3,size(g,3),size(g,4));
% cdgx(:,:,:,1:size(g,3)-1,1:size(g,4)-1)=dgx;

%% Interpolation of g, gi (g inverse)  and dgx (derivative of g along x (x=x1,x2,x3))
% [X,Y]=meshgrid(1:size(g,4),1:size(g,3));
% [Xq,Yq] = meshgrid(1:0.25:size(g,4), 1:0.25:size(g,3));
% for ei=1:3
%     for ej=1:3
%         g_inter(ei,ej,:,:)=interp2(X,Y,squeeze(g(ei,ej,:,:)),Xq,Yq);
%         gi_inter(ei,ej,:,:)=interp2(X,Y,squeeze(gi(ei,ej,:,:)),Xq,Yq);
%         dgx_inter(ei,ej,1,:,:)=interp2(X,Y,squeeze(cdgx(ei,ej,1,:,:)),Xq,Yq);
%         dgx_inter(ei,ej,2,:,:)=interp2(X,Y,squeeze(cdgx(ei,ej,2,:,:)),Xq,Yq);
%         dgx_inter(ei,ej,3,:,:)=interp2(X,Y,squeeze(cdgx(ei,ej,3,:,:)),Xq,Yq);
%     end
% end

%% Christoffel Symbol
lopsz=3;
for x=2:size(gi_inter,3)-1
    for y=2:size(gi_inter,4)-1
        for i=1:lopsz
            for j=1:lopsz
                for k=1:lopsz
                    c=0;
                    for l=1:lopsz
                        c = c + (1/2)*gi_inter(k,l,x,y)*(dgx_inter(i,l,j,x,y)+dgx_inter(j,l,i,x,y)-dgx_inter(i,j,l,x,y));
                    end
                    crst_inter(i,j,k,x,y)=c;
                end
            end
        end
    end
end

%% Central derivative of Christoffel Symbol along X and Y
for x=3:size(gi_inter,3)-2
    for y=3:size(gi_inter,4)-2
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
% dcrstx_inter=zeros(3,3,3,3,size(gi_inter,3),size(gi_inter,4));
% dcrstx_inter(:,:,:,:,1:size(gi_inter,3)-1,1:size(gi_inter,4)-1)=temp_dcrstx_inter;
dcrstx_inter=temp_dcrstx_inter;
%% Ricci Tensor
for x=3:size(gi_inter,3)-2
    for y=3:size(gi_inter,4)-2
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
RS_by_RicTensor= (sum(sum(gi_inter(:,:,1:size(gi_inter,3)-2,1:size(gi_inter,3)-2).*RicTensor_inter))); % compute Ricci Scalar from Ricci Tensor
 RS(1:size(gi_inter,3)-2,1:size(gi_inter,4)-2)=RS_by_RicTensor;
figure; surf(RS(3:size(gi_inter,3)-2,3:size(gi_inter,4)-2));
% figure; imagesc(RS);
% RS_All(1:size(gi_inter,3),1:size(gi_inter,4),ii)=RS_by_RicTensor;
% end
% RS_Sum=zeros(size(gi_inter,3),size(gi_inter,4));
% for ii=1:max_num
%     RS_Sum=RS_Sum+ squeeze(RS_All(:,:,ii));
% end
% figure; surf(RS_Sum);