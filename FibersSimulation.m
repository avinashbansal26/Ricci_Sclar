function [Tensor_ODF]=FibersSimulation()

% global T9_2D;
order=4;%In this demo we compute 4th-order tensorsfor ii = 1:9
%     
    S0=1.0;
UnitVectors;
n_grad=81;
GradientOrientations=[1 0 0; g(1:n_grad,:)];
b_value=[10;ones(n_grad,1)*1500];
lsize=32;
bvalue=1500;
S=ones(lsize,lsize,1,size(GradientOrientations,1));
 f1count=0;
 f2count=0;
 finter_count=0;
 fbckg_count=0;

%  rr=randi([1 20],10,1);

for x=1:lsize
    for y=1:lsize
        
        for i=2:size(GradientOrientations,1)
             
          f1=0;f2=0;f3=0;  flag=0;
          
           %%Shape 1 Intersection Horizontal and Vertical
           if y<26 & y>22 & x<32 & x>1
               fiber2=[ 1 0 0];
                 f2=1;
           end
           
           
           %Define the underlying fiber orientation of fiber bundle #2 (i.e. we have 2-fiber crossing)
%            if x<19 & x>14 & y<32 & y>1
%                fiber2=[ 0 1 0];
%                f2=1;
%            end
%      Shape 2, 3, 4: Intersecting fibers circle, line having upto 7 regions     
%% geodesic1
% if(x<29 & y>3 & y<29)
%     xx=8-(x/4);
%     n=.35;yy=n*y;%0.5,1,n=1/4 0.35 each!!,0.25 this one also
%     s=1;
%     if((xx-sin(yy))>2 & (xx-sin(yy))<4)
%         
%         v=[-(1/s)*cos(yy) xx 0];%v=v/sqrt(v*v');
%         
%         fiber2=v;    f2=1;
%         %
%     end
% end
%% geodesic 2 U shape
if y<32 & y>1 & x<32 & x>1
if( y>1 & y<32) 
yy=(y-3);
xx=(x-16);

% 
if xx*xx+yy*yy>4*4 & xx*xx+yy*yy<7*7
    theta=atan(yy/xx);
    v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
    fiber2=v;    f2=1;
end

end
end
% % % and a line..to create  deep rev-U
%  if y>1 & y<2 & x>9 & x<13
%      fiber2=[ 0 1 0];
%      f2=1;
%  end
%  if y>1 & y<2 & x>19 & x<23
%      fiber2=[ 0 1 0];
%      f2=1;
%  end
%% temp code you can delete this code intersection of two lines
% 
if y<32 & y>1 & x<32 & x>1
    theta=45*pi/180; % pair (0.8,0.6) and (0.8 amd 0.7) 0.6 This one also
            m=tan(theta);
            if y-m*x>0 & y-m*x<4  % m=0.3,0.7,0.4--0.1,0.2 fails
                v=[ cos(theta) sin(theta) 0 ];%v=v/sqrt(v*v');
                fiber1=v;
                f1=1;               
            end
            
            theta=45*pi/180; % pair (0.8,0.6) and (0.8 amd 0.7) 0.6 This one also
            m=tan(theta);
            if y-m*x>11 & y-m*x<16  % m=0.3,0.7,0.4--0.1,0.2 fails
                v=[ cos(theta) sin(theta) 0 ];%v=v/sqrt(v*v');
                fiber1=v;
                f1=1;               
            end
end
            if y<32 & y>1 & x<32 & x>1
            theta=135*pi/180; % pair (0.8,0.6) and (0.8 amd 0.7) 0.6 This one also
            m=tan(theta);
           if y-m*x>32 & y-m*x<37  % m=0.3,0.7,0.4--0.1,0.2 fails           
           v=[ cos(theta) sin(theta) 0 ];
           fiber3=v;
           f3=1;
            end
            end
% 

%% geodesic 3 Cubic not working

%         xx=((16-x)^(-2));
%         yy=abs(16-y);
%         yy-xx^3
%         if((xx^3-yy)>2 && (xx^3-yy)<16)
%            theta=atan((yy)/(xx));
%            v=[cos(theta)  sin(theta) 0];v=v/sqrt(v*v'); 
% %             v=[ 3*(xx)^2   yy 0];v=v/sqrt(v*v');
%             
%             fiber2=v;    f2=1;
%         end

%% and a line..
% if y<12 & y>8 %& x>19 & x<10
%     fiber2=[ 1 0 0];
%     f2=1;
% end
%% S-shape Complex Curved Geodesic 4 for testing Adjugate sigmoid etc..

% if(x<16)
%     xx=(16-x);
%     if( y>7 & y<21)
%         yy=(y-7);
%         
%         
%         if xx*xx+yy*yy>9*9 & xx*xx+yy*yy<14*14
%             theta=atan(yy/xx);
%             v=[sin(theta)  cos(theta) 0];v=v/sqrt(v*v');
%             fiber2=v;    f2=1;
%         end
%     end
% end
% if(x>=16 & x<32)
%     xx=(16-x);
%     if( y>12 & y<31)
%         yy=29-y;
%         
%         if xx*xx+yy*yy>9*9 & xx*xx+yy*yy<14*14
%             theta=atan(yy/xx);
%             v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
%             fiber2=v;    f2=1;
%         end
%     end
% end
% if y>2 & y<8 & x>2 & x<8
%     fiber2=[ 0 1 0];
%     f2=1;
% end
%% geodesic 5 to create more smooth flow cubic like shape
% % 
% if(x<=13)
% xx=(13-x);
% if(y>2 & y<31) 
% yy=11-y;
% if xx*xx+yy*yy>5*5 & xx*xx+yy*yy<10*10
%     theta=atan(yy/xx);
%     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
%     fiber2=v;    f2=1;
% end
% end 
% end 
% %%%
% if(x>=22)
% xx=(22-x);
% if(y>=16 & y<31) 
% yy=23-y;
% if xx*xx+yy*yy>3*3 & xx*xx+yy*yy<8*8
%     theta=atan(yy/xx);
%     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
%     fiber2=v;    f2=1;
% end
% end 
% end
% if y>15 & y<21 & x>11 & x<22
%      fiber2=[ 1 0 0]; 
%      f2=1;
%  end

%% %%%%%%%%%%%%%%%%%%%%%%
%            
        %%      Helix Shape
%          yy=abs(9-y);
%          xx=abs(9-x);
%         
%         if((yy-xx*tan(sqrt(xx*xx+yy*yy)))>2 & (yy-xx*tan(sqrt(xx*xx+yy*yy)))<16)
%             theta=atan((9-y)/(9-x));
%             v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
%             fiber1=v;    f1=1;
%         end   
        %%        Sine Shape  It works
%         if(x>5 && x<25)
%         yy=6-(y/5);%yy=8-(y/1.5);
%         n=.45;xx=n*x;%0.5,1,n=1/4
%         s=1;
%         if((yy-sin(xx/s))>3 & (yy-sin(xx/s))<4.4)
%             
%             v=[yy -(1/s)*cos(xx) 0];v=v/sqrt(v*v');
%             fiber1=v;    f1=1;
%             
%         end
%       end   
%         xx=8-(x/2.7);
%         n=.25;yy=n*y;%0.5,1,n=1/4 0.35 each!!,0.25 this one also
%         s=1;
%         if((xx-sin(yy))>2 & (xx-sin(yy))<3.5)
%             
%             v=[-(1/s)*cos(yy) xx 0];%v=v/sqrt(v*v');
%             
%             fiber2=v;    f2=1;
% %             
%         end
%         %%       parabola works   
%         yy=(y);
%         xx=x-8;
%         
%         if((yy-(xx^2))>0.5 & (yy-(xx^2))<10)
%            theta=atan((yy)/(xx));
% %            v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v'); 
%             v=[yy 2*(xx) 0];v=v/sqrt(v*v');
%             
%             fiber2=v;    f2=1;
%         end
        



%  %%        Sine Shape  It works
% %                    if(x>5 && x<25)
%                         yy=7-(y/5);%yy=8-(y/1.5);
%                         n=.1;xx=n*x;%0.5,1,n=1/4
%                         s=1;
%                         if((yy-sin(xx/s))>3 & (yy-sin(xx/s))<4.4)
%                 
%                             v=[yy -(1/s)*cos(xx) 0];v=v/sqrt(v*v');
%                             fiber1=v;    f1=1;
%                 
%                         end
% %                       end
%                         xx=6-(x/5);
%                         n=.1;yy=n*y;%0.5,1,n=1/4 0.35 each!!,0.25 this one also
%                         s=1;
%                         if((xx-sin(yy))>2 & (xx-sin(yy))<3.5)
%                 
%                             v=[-(1/s)*cos(yy) xx 0];%v=v/sqrt(v*v');
%                 
%                             fiber2=v;    f2=1;
%                 %
%                         end
          %% Shape 2, 3, 4: Intersecting fibers circle, line having upto 7 regions
%             if x*x+y*y>10*10 & x*x+y*y<14*14% 25*25,32*32
% %                 theta=atan(y/x);s=1;
% %                 theta=atan(tan(theta)/s);
%                 v=[y/x -1 0];v=v/sqrt(v*v');
% %                  v=[sin(theta) -cos(theta) 0];v=v/sqrt(v*v');
%                  
%                 fiber2=v;
%                 f2=1;
%             end


%             m=-0.2;%0.8 0.2 for one result sine and line 0.6 is merging
%             if y-m*x>20 & y-m*x<25  % m=0.3,0.7,0.4--0.1,0.2 fails
%                 theta=atan(m); 
%                 v=[  cos(theta) sin(theta) 0 ];%v=v/sqrt(v*v');
%                 fiber2=v;
%                 f2=1;               
%             end
%             m=0.3; % pair (0.8,0.6) and (0.8 amd 0.7) 0.6 This one also
%            if y-m*x>10 & y-m*x<15  % m=0.3,0.7,0.4--0.1,0.2 fails
%            theta=atan(m); 
%            v=[  cos(theta)  sin(theta) 0 ];
%            fiber1=v;
%            f1=1;
%            end




%             if ((12-x)*(12-x)+y*y)>12*12 & ((10-x)*(10-x)+y*y)<15*15% 7, 11 fails
%                      
%                 v=[ y (12-x) 0 ];v=v/sqrt(v*v');%22-x
%                 fiber2=v;
%                 f2=1;
%             end   
    %% Note 0.8 vs 0.6,0.5,0.4 results are important
%               m=0.01;%0.6 0.2 for one result sine and line 0.6 is merging
%             if y-m*x>10 & y-m*x<14  % m=0.3,0.7,0.4--0.1,0.2 fails
%                 theta=atan(m); s=1;
% %                 theta=atan(tan(theta)/s);
%                  v=[ cos(theta) sin(theta) 0 ];%v=v/sqrt(v*v');
% %                 thet=0*pi/180;
% %                 thet=0;
% %                 v1=[ cos(thet)   -sin(thet)      0;
% %                      sin(thet)   cos(thet)     0;
% %                      0   0  1;]*v';
% %                  v1=v1'/sqrt(v1'*v1);
%                 fiber1=v;
%                 f1=1;
%             end
%             if y>2 & y<5  % 5 10
%                  v=[ 1 0 0 ];
%                 fiber3=v;
%                 f3=1;
%             end
%             flag=0;


            %%Intersections
            if(f1==1 & f2==0 & f3==0)
                pf1=1;
                pf2=0;
                pf3=0;
%                 r = rand(1, 2); % Start with 3 random numbers that don't sum to 1.
%                  r = r / sum(r);  % Normalize so the sum is 1.
%                  pf1=r(1);
%                  pf2=0;
%                  pf3=0;
                f1count=f1count+1;
            elseif(f1==0 & f2==1 & f3==0)   
                pf1=0;
                pf2=1;
                pf3=0;
%                 r = rand(1, 2); % Start with 3 random numbers that don't sum to 1.
%                  r = r / sum(r);  % Normalize so the sum is 1.
%                  pf1=0;
%                  pf2=r(1);
%                  pf3=0;
                f2count=f2count+1;
            elseif(f1==0 & f2==0 & f3==1)   
                pf1=0;
                pf2=0;
                pf3=1;
               
            elseif(f1==1 & f2==1 & f3==0)   
                
                pf1=0.5;
                pf2=0.5;
                pf3=0;
%                  r = rand(1, 2); % Start with 3 random numbers that don't sum to 1.
%                  r = r / sum(r);  % Normalize so the sum is 1.
%                  pf1=r(1);
%                  pf2=r(2);
%                  pf3=0;
               finter_count=finter_count+1;
            elseif(f1==0 & f2==1 & f3==1)   
                pf1=0;
                pf2=0.5;
                pf3=0.5;
%                  r = rand(1, 2); % Start with 3 random numbers that don't sum to 1.
%                  r = r / sum(r);  % Normalize so the sum is 1.
%                  pf1=r(1);
%                  pf2=r(2);
%                  pf3=0;
               
            elseif(f1==1 & f2==0 & f3==1)   
                pf1=0.5;
                pf2=0;
                pf3=0.5;
              
            elseif(f1==1 & f2==1 & f3==1) 
                pf1=0.333333;
                pf2=0.333333;
                pf3=0.333333;   
            else   
                pf1=0;
                pf2=0;
                pf3=0;
%                  r = rand(1, 3); % Start with 3 random numbers that don't sum to 1.
%                  r = r / sum(r);  % Normalize so the sum is 1.
%                  pf1=r(1);
%                  pf2=r(2);
%                  pf3=r(3);
                
               fiber1=[1 0 0];
               fiber2=[0 1 0]; 
               fiber3=[0 0 1]; 
               flag=1;
               fbckg_count=fbckg_count+1;
%                 pf1=0.8;%0.8,0.1,0.1 circle
%                 pf2=0.1;
%                 pf3=0.1;
%                fiber1=[1 0 0];
%                fiber2=[0 1 0]; 
%                fiber3=[0 0 1]; 
               
            end
                
                              
                              
             
           if(flag==0)  
           S(x,y,1,i)=S(x,y,1,1)*(pf1*SimulateDWMRI(fiber1,GradientOrientations(i,:))+ ...
                      +pf2*SimulateDWMRI(fiber2,GradientOrientations(i,:))+...
                      +pf3*SimulateDWMRI(fiber3,GradientOrientations(i,:)));
           else
               S(x,y,1,i)=(0.01);
           end
           
%             S(x,y,1,i)=S(x,y,1,1)*(pf1*SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+ pf2*SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:)));
          
        end
    end    
 end
f1count=f1count/81;
f2count=f2count/81;
finter_count=finter_count/81;
fbckg_count=fbckg_count/81;
fregion=[f1count,f2count,finter_count,fbckg_count];
%% Adding Riccian Noise
% sigma=0.00;%0.02[0.01-0.09]in tha paper
% for x=1:lsize
%     for y=1:lsize
%         for i=1:size(GradientOrientations,1)
%             S(x,y,1,i)=sqrt((S(x,y,1,i)+sigma*randn(1))^2+(sigma*randn(1))^2);
%         end
%     end
% end


%Construct all possible monomials of a specific order
G=constructMatrixOfMonomials(GradientOrientations, order); %computes G from section 5.1 (ISBI'10)
%Construct set of polynomial coefficients C
C=constructSetOf321Polynomials(order)'; %computes C from section 5.1 (ISBI'10)
P=G*C;
P=[-diag(b_value)*P ones(size(GradientOrientations,1),1)];
%%For CT ODF
BG=constructMatrixOfIntegrals(GradientOrientations, order, 100);
B=BG*C;

f=sqrt(2);
for i=1:size(S,1)
    
    for j=1:size(S,2)
       
        y=squeeze(log(S(i,j,1,:)));
        x=lsqnonneg(P, y);
        UniqueTensorCoefficients = C * x([1:321]);
        
        TensorCoefficients(:,i,j)=UniqueTensorCoefficients;
         D=UniqueTensorCoefficients;
         %% For CT ODF
         x1=lsqnonneg(B, exp(-bvalue*G*D));
         CT_Coeff = C * x1;
         Tensor_ODF(:,i,j)=CT_Coeff;
         TD=CT_Coeff;
         %% CT ODF ends
        
%         T6D(i,j,:,:)=[A1111,   A1122,   A1133,   f*A1112,  f*A1113,  f*A1123;
%                       A1122,   A2222,   A2233,   f*A2212,  f*A2213,  f*A2223;
%                       A1133,   A2233,   A3333,   f*A3312,  f*A3313,  f*A3323;
%                       f*A1112, f*A2212, f*A3312, 2*A1122,  2*A1123,  2*A2213;
%                       f*A1113, f*A2213, f*A3313, 2*A1123,  2*A1133,  2*A3312;
%                       f*A1123, f*A2223, f*A3323, 2*A2213,  2*A3312,  2*A2233;];
%                 chol(squeeze(T6D(i,j,:,:)))
        %
        dxx=(3/35)*(9*D(15)+8*(1/6)*D(12)+ 8*(1/6)*D(10)-D(5)-D(1)-(1/6)*D(3));
        dyy=(3/35)*(9*D(5)+ 8*(1/6)*D(12)+8*(1/6)*D(3)-D(15)-D(1)-(1/6)*D(10));
        dzz=(3/35)*(9*D(1)+8*(1/6)*D(10)+ 8*(1/6)*D(3)-D(15)-D(5)-(1/6)*D(12));
        dxy=(6/7)*((1/4)*D(14)+(1/4)*D(9)+(1/12)*D(7));
        dxz=(6/7)*((1/4)*D(13)+(1/4)*D(9)+(1/12)*D(7));
        dyz=(6/7)*((1/4)*D(4)+(1/4)*D(2)+(1/12)*D(11));
        % dxx=(3/35)*(9*D(15)+8*D(12)+ 8*D(10)-D(5)-D(1)-2*D(3));
        % dyy=(3/35)*(9*D(5)+ 8*D(12)+8*D(3)-D(15)-D(1)-2*D(10));
        % dzz=(3/35)*(9*D(1)+ 8*D(10)+ 8*D(3)-D(15)-D(5)-2*D(12));
        % dxy=(6/7)*(D(14)+D(9)+D(7));
        % dxz=(6/7)*(D(13)+D(9)+D(7));
        % dyz=(6/7)*(D(4)+D(2)+(11));
        
        %%From Ozarslan paper
        T6_2D(i,j,:,:)=[dxx,0.5*dxy,0.5*dxz;
            0.5*dxy,dyy,0.5*dyz;
            0.5*dxz,0.5*dyz,dzz;];
        tLowDim(:,:,i,j)=T6_2D(i,j,:,:);
        % chol(squeeze(T6_2D(i,j,:,:)));
        T6D(i,j,:,:)= [D(15),         (1/6)*D(12),    (1/6)*D(10),   (1/4)*f*D(14), (1/4)*f*D(13),  (1/12)*f*D(11);
                      (1/6)*D(12),     D(5),          (1/6)*D(3),    (1/4)*f*D(9),  (1/12)*f*D(8),  (1/4)*f*D(4);
                      (1/6)*D(10),    (1/6)*D(3),     D(1),          (1/12)*f*D(7), (1/4)*f*D(6),   (1/4)*f*D(2);
                      (1/4)*f*D(14),  (1/4)*f*D(9),   (1/12)*f*D(7), (1/6)*2*D(12), (1/12)*2*D(11), (1/12)*2*D(8);
                      (1/4)*f*D(13),  (1/12)*f*D(8),  (1/4)*f*D(6),  (1/12)*2*D(11),(1/6)*2*D(10),  (1/12)*2*D(7);
                      (1/12)*f*D(11), (1/4)*f*D(4),   (1/4)*f*D(2),  (1/12)*2*D(8), (1/12)*2*D(7),  (1/6)*2*D(3);];
        
        T9D(i,j,:,:)=[  D(15),  (1/4)*D(14),  (1/4)*D(13),  (1/4)*D(14),  (1/6)*D(12), (1/12)*D(11),  (1/4)*D(13),   (1/12)*D(11), (1/6)*D(10);
                        (1/4)*D(14),  (1/6)*D(12),  (1/12)*D(11), (1/6)*D(12),  (1/4)*D(9),  (1/12)*D(8),   (1/12)* D(11), (1/12)*D(8),  (1/12)*D(7);
                        (1/4)*D(13),  (1/12)*D(11), (1/6)*D(10),  (1/12)*D(11), (1/12)*D(8), (1/12)*D(7),   (1/6)*D(10),   (1/12)*D(7),  (1/4)*D(6);
            
                        (1/4)*D(14),  (1/6)*D(12),  (1/12)*D(11), (1/6)*D(12),  (1/4)*D(9),   (1/12)*D(8),  (1/12)*D(11),  (1/12)*D(8),  (1/12)*D(7);
                        (1/6)*D(12),  (1/4)*D(9),   (1/12)*D(8),  (1/4)*D(9),    D(5),        (1/4)*D(4),   (1/12)*D(8),   (1/4)*D(4),   (1/6)*D(3);
                        (1/12)*D(11), (1/12)*D(8),  (1/12)*D(7),  (1/12)*D(8),   (1/4)*D(4),  (1/6)*D(3),   (1/12)*D(7),   (1/6)*D(3),   (1/4)*D(2);
            
                        (1/4)*D(13),  (1/12)*D(11), (1/6)*D(10),  (1/12)*D(11),  (1/12)*D(8), (1/12)*D(7),  (1/6)*D(10),   (1/12)*D(7),  (1/4)*D(6);
                        (1/12)* D(11),(1/12)*D(8),  (1/12)*D(7),  (1/12)*D(8),   (1/4)*D(4),  (1/6)*D(3),   (1/12)*D(7),   (1/6)*D(3),   (1/4)*D(2);
                        (1/6)*D(10),  (1/12)*D(7),  (1/4)*D(6),   (1/12)*D(7),   (1/6)*D(3),  (1/4)*D(2),   (1/4)*D(6),    (1/4)*D(2),   D(1);];
        T9_CT(i,j,:,:)=[  TD(15),  (1/4)*TD(14),  (1/4)*TD(13),  (1/4)*TD(14),  (1/6)*TD(12), (1/12)*TD(11),  (1/4)*TD(13),   (1/12)*TD(11), (1/6)*TD(10);
                        (1/4)*TD(14),  (1/6)*TD(12),  (1/12)*TD(11), (1/6)*TD(12),  (1/4)*TD(9),  (1/12)*TD(8),   (1/12)* TD(11), (1/12)*TD(8),  (1/12)*TD(7);
                        (1/4)*TD(13),  (1/12)*TD(11), (1/6)*TD(10),  (1/12)*TD(11), (1/12)*TD(8), (1/12)*TD(7),   (1/6)*TD(10),   (1/12)*TD(7),  (1/4)*TD(6);
            
                        (1/4)*TD(14),  (1/6)*TD(12),  (1/12)*TD(11), (1/6)*TD(12),  (1/4)*TD(9),   (1/12)*TD(8),  (1/12)*TD(11),  (1/12)*TD(8),  (1/12)*TD(7);
                        (1/6)*TD(12),  (1/4)*TD(9),   (1/12)*TD(8),  (1/4)*TD(9),    TD(5),        (1/4)*TD(4),   (1/12)*TD(8),   (1/4)*TD(4),   (1/6)*TD(3);
                        (1/12)*TD(11), (1/12)*TD(8),  (1/12)*TD(7),  (1/12)*TD(8),   (1/4)*TD(4),  (1/6)*TD(3),   (1/12)*TD(7),   (1/6)*TD(3),   (1/4)*TD(2);
            
                        (1/4)*TD(13),  (1/12)*TD(11), (1/6)*TD(10),  (1/12)*TD(11),  (1/12)*TD(8), (1/12)*TD(7),  (1/6)*TD(10),   (1/12)*TD(7),  (1/4)*TD(6);
                        (1/12)* TD(11),(1/12)*TD(8),  (1/12)*TD(7),  (1/12)*TD(8),   (1/4)*TD(4),  (1/6)*TD(3),   (1/12)*TD(7),   (1/6)*TD(3),   (1/4)*TD(2);
                        (1/6)*TD(10),  (1/12)*TD(7),  (1/4)*TD(6),   (1/12)*TD(7),   (1/6)*TD(3),  (1/4)*TD(2),   (1/4)*TD(6),    (1/4)*TD(2),   TD(1);];
        T6_CT(i,j,:,:)= [TD(15),         (1/6)*TD(12),    (1/6)*TD(10),   (1/4)*f*TD(14), (1/4)*f*TD(13),  (1/12)*f*TD(11);
                      (1/6)*TD(12),     TD(5),          (1/6)*TD(3),    (1/4)*f*TD(9),  (1/12)*f*TD(8),  (1/4)*f*TD(4);
                      (1/6)*TD(10),    (1/6)*TD(3),     TD(1),          (1/12)*f*TD(7), (1/4)*f*TD(6),   (1/4)*f*TD(2);
                      (1/4)*f*TD(14),  (1/4)*f*TD(9),   (1/12)*f*TD(7), (1/6)*2*TD(12), (1/12)*2*TD(11), (1/12)*2*TD(8);
                      (1/4)*f*TD(13),  (1/12)*f*TD(8),  (1/4)*f*TD(6),  (1/12)*2*TD(11),(1/6)*2*TD(10),  (1/12)*2*TD(7);
                      (1/12)*f*TD(11), (1/4)*f*TD(4),   (1/4)*f*TD(2),  (1/12)*2*TD(8), (1/12)*2*TD(7),  (1/6)*2*TD(3);];
        
%          T6D_CT(i,j,:,:)= [  D_CT(15),         (1/6)*D_CT(12),    (1/6)*D_CT(10),   (1/4)*f*D_CT(14), (1/4)*f*D_CT(13),  (1/12)*f*D_CT(11);
%                             (1/6)*D_CT(12),     D_CT(5),          (1/6)*D_CT(3),    (1/4)*f*D_CT(9),  (1/12)*f*D_CT(8),  (1/4)*f*D_CT(4);
%                             (1/6)*D_CT(10),    (1/6)*D_CT(3),     D_CT(1),          (1/12)*f*D_CT(7), (1/4)*f*D_CT(6),   (1/4)*f*D_CT(2);
%                             (1/4)*f*D_CT(14),  (1/4)*f*D_CT(9),   (1/12)*f*D_CT(7), (1/6)*2*D_CT(12), (1/12)*2*D_CT(11), (1/12)*2*D_CT(8);
%                             (1/4)*f*D_CT(13),  (1/12)*f*D_CT(8),  (1/4)*f*D_CT(6),  (1/12)*2*D_CT(11),(1/6)*2*D_CT(10),  (1/12)*2*D_CT(7);
%                             (1/12)*f*D_CT(11), (1/4)*f*D_CT(4),   (1/4)*f*D_CT(2),  (1/12)*2*D_CT(8), (1/12)*2*D_CT(7),  (1/6)*2*D_CT(3);];
%         %    T9fromFourthOrder(i,j,:,:)=[squeeze(W(i,j,:,:,1,1), squeeze(W(i,j,:,:,1,2)), squeeze(W(i,j,:,:,1,3));
%         %                                squeeze(W(i,j,:,:,2,1), squeeze(W(i,j,:,:,2,2) , squeeze(W(i,j,:,:,2,3));
%         %                                squeeze(W(i,j,:,:,3,1), squeeze(W(i,j,:,:,3,2) , squeeze(W(i,j,:,:,3,3))];
%          T9D_CT(i,j,:,:)=[        D_CT(15),  (1/4)*D_CT(14),  (1/4)*D_CT(13),  (1/4)*D_CT(14),  (1/6)*D_CT(12), (1/12)*D_CT(11),  (1/4)*D_CT(13),   (1/12)*D_CT(11), (1/6)*D_CT(10);
%                                  (1/4)*D_CT(14),  (1/6)*D_CT(12),  (1/12)*D_CT(11), (1/6)*D_CT(12),  (1/4)*D_CT(9),  (1/12)*D_CT(8),   (1/12)* D_CT(11), (1/12)*D_CT(8),  (1/12)*D_CT(7);
%                                  (1/4)*D_CT(13),  (1/12)*D_CT(11), (1/6)*D_CT(10),  (1/12)*D_CT(11), (1/12)*D_CT(8), (1/12)*D_CT(7),   (1/6)*D_CT(10),   (1/12)*D_CT(7),  (1/4)*D_CT(6);
%             
%                                  (1/4)*D_CT(14),  (1/6)*D_CT(12),  (1/12)*D_CT(11), (1/6)*D_CT(12),  (1/4)*D_CT(9),   (1/12)*D_CT(8),  (1/12)*D_CT(11),  (1/12)*D_CT(8),  (1/12)*D_CT(7);
%                                  (1/6)*D_CT(12),  (1/4)*D_CT(9),   (1/12)*D_CT(8),  (1/4)*D_CT(9),    D_CT(5),        (1/4)*D_CT(4),   (1/12)*D_CT(8),   (1/4)*D_CT(4),   (1/6)*D_CT(3);
%                                  (1/12)*D_CT(11), (1/12)*D_CT(8),  (1/12)*D_CT(7),  (1/12)*D_CT(8),   (1/4)*D_CT(4),  (1/6)*D_CT(3),   (1/12)*D_CT(7),   (1/6)*D_CT(3),   (1/4)*D_CT(2);
%             
%                                  (1/4)*D_CT(13),  (1/12)*D_CT(11), (1/6)*D_CT(10),  (1/12)*D_CT(11),  (1/12)*D_CT(8), (1/12)*D_CT(7),  (1/6)*D_CT(10),   (1/12)*D_CT(7),  (1/4)*D_CT(6);
%                                  (1/12)* D_CT(11),(1/12)*D_CT(8),  (1/12)*D_CT(7),  (1/12)*D_CT(8),   (1/4)*D_CT(4),  (1/6)*D_CT(3),   (1/12)*D_CT(7),   (1/6)*D_CT(3),   (1/4)*D_CT(2);
%                                  (1/6)*D_CT(10),  (1/12)*D_CT(7),  (1/4)*D_CT(6),   (1/12)*D_CT(7),   (1/6)*D_CT(3),  (1/4)*D_CT(2),   (1/4)*D_CT(6),    (1/4)*D_CT(2),   D_CT(1);];
        
                  chol(squeeze(T6D(i,j,:,:)));
     
%          Tensor_ODF(:,i,j)=D;
    end
end

%% Till here We got the 6x6 tensors
%%
figure;
plotTensors(TensorCoefficients,1,[321  1]);
% figure;plotTensors(TensorCoefficients(:,16,16),1,[321 1] );
% figure;
% plotTensors(tensor_CT,1,[321 1]);

figure;plotTensors(Tensor_ODF,1,[321  1]);


figure;plotDTI(tLowDim,.005);
end
