        function [Tensor_ODF]=geodesic_simulation_cross()

    % global T9_2D;
    order=4;%In this demo we compute 4th-order tensorsfor ii = 1:9
    %
    S0=1.0;
    UnitVectors;
    n_grad=81;
    GradientOrientations=[1 0 0;g([1:n_grad],:)];
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

                f1=0;f2=0;f3=0;

                %% Shape 1 Intersection Horizontal and Vertical
%                            if  ((x>=1 && x<=32) && (y<19 && y>13))
%                                fiber1=[ 1 0 0];
%                                  f1=1;
%                            end                
%                 %            %Define the underlying fiber orientation of fiber bundle #2 (i.e. we have 2-fiber crossing)
%                            if ((x<19 && x>12) && (y>12 && y<31))
%                                fiber2=[ 0 1 0];
%                                f2=1;
%                            end
                          
                %%      Shape 2, 3, 4: Intersecting fibers circle, line having upto 7 regions
                %%  geodesic1

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
 %% two mering fibers in C shape              
%                 if( y>3 & y<19)
%                 yy=(y-5);
%                 xx=(x-16);
%                 if xx*xx+yy*yy>8*8 & xx*xx+yy*yy<13*13
%                     theta=atan(yy/xx);
%                     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
%                     fiber1=v;    f1=1;
%                 end
%                 end
%                 if( y>16 & y<30)
%                 yy=(y-27);
%                 xx=(x-16);
%                 if xx*xx+yy*yy>8*8 & xx*xx+yy*yy<13*13
%                     theta=atan(yy/xx);
%                     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
%                     fiber2=v;    f2=1;
%                 end
%                 end
%                 
                
                

                %%   geodesic 2 U shape
                
                if( y>1 & y<32)
                yy=(y-16);
                xx=(x-16);
                if xx*xx+yy*yy>11*11 & xx*xx+yy*yy<15*15
                    theta=atan(yy/xx);
                    v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                    fiber2=v;    f2=1;
                end
                end

                % % % and a line..to create  deep rev-U
%                  if x>5 & x<28 & y>17 & y<23
%                      fiber2=[ 1 0 0];
%                      f2=1;
%                  end
%                  if y>2 & y<19 & x>23 & x<29
%                      fiber2=[ 0 1 0];
%                      f2=1;
%                  end
                

%%  temp Code you can use this code

            theta=55*pi/180; % pair (0.8,0.6) and (0.8 amd 0.7) 0.6 This one also
            m=tan(theta);
            if y-m*x>0 & y-m*x<7  % m=0.3,0.7,0.4--0.1,0.2 fails
                v=[ cos(theta) sin(theta) 0 ];%v=v/sqrt(v*v');
                fiber1=v;
                f1=1;               
            end
            theta=118*pi/180; % pair (0.8,0.6) and (0.8 amd 0.7) 0.6 This one also
            m=tan(theta);
           if y-m*x>54 & y-m*x<62  % m=0.3,0.7,0.4--0.1,0.2 fails
           
           v=[ cos(theta) sin(theta) 0 ];
           fiber3=v;
           f3=1;
           end

%%  Complex spiral Shape
                % if( y>=18 & y<32)
                % yy=(y-18);
                % xx=(x-19);
                % if xx*xx+yy*yy>2*2 & xx*xx+yy*yy<5*5
                %     theta=atan(yy/xx);
                %     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                %     fiber2=v;    f2=1;
                % end
                % end
                %
                % if( y>1 & y<=18)
                % yy=(y-18);
                % xx=(x-16);
                % if xx*xx+yy*yy>5*5 & xx*xx+yy*yy<8*8
                %     theta=atan(yy/xx);
                %     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                %     fiber2=v;    f2=1;
                % end
                % end
                %
                % if( y>=18 & y<32)
                % yy=(y-18);
                % xx=(x-19);
                % if xx*xx+yy*yy>8*8 & xx*xx+yy*yy<11*11
                %     theta=atan(yy/xx);
                %     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                %     fiber2=v;    f2=1;
                % end
                % end
                %
                % if( y>1 & y<=18)
                % yy=(y-18);
                % xx=(x-16);
                % if xx*xx+yy*yy>11*11 & xx*xx+yy*yy<14*14
                %     theta=atan(yy/xx);
                %     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                %     fiber2=v;    f2=1;
                % end
                % end

                %%   3's Shape (Easy Three)
%                 if( x>=16 & x<32)
%                     yy=(y-12);
%                     xx=(x-14);
%                     if xx*xx+yy*yy>4*4 & xx*xx+yy*yy<7*7
%                         theta=atan(yy/xx);
%                         v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
%                         fiber2=v;    f2=1;
%                     end
%                 end
% 
%                 if( x>=16 & x<32)
%                     yy=(y-22);
%                     xx=(x-14);
%                     if xx*xx+yy*yy>4*4 & xx*xx+yy*yy<7*7
%                         theta=atan(yy/xx);
%                         v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
%                         fiber1=v;    f1=1;
%                     end
%                 end


                %%   3's Shape (complex Three)
                % if( x>=13 & x<32)
                % yy=(y-12);
                % xx=(x-14);
                % if xx*xx+yy*yy>4*4 & xx*xx+yy*yy<7*7
                %     theta=atan(yy/xx);
                %     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                %     fiber2=v;    f2=1;
                % end
                % end
                %
                % if( x>=13 & x<32)
                % yy=(y-22);
                % xx=(x-14);
                % if xx*xx+yy*yy>4*4 & xx*xx+yy*yy<7*7
                %     theta=atan(yy/xx);
                %     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                %     fiber1=v;    f1=1;
                % end
                % end

                %%   Fork 3's Shape (Three)
                % if( x>=9 & x<32)
                % yy=(y-12);
                % xx=(x-10);
                % if xx*xx+yy*yy>4*4 & xx*xx+yy*yy<7*7
                %     theta=atan(yy/xx);
                %     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                %     fiber2=v;    f2=1;
                % end
                % end
                %
                % if( x>=9 & x<32)
                % yy=(y-22);
                % xx=(x-10);
                % if xx*xx+yy*yy>4*4 & xx*xx+yy*yy<7*7
                %     theta=atan(yy/xx);
                %     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                %     fiber1=v;    f1=1;
                % end
                % end
                %
                % if y<19 & y>15 & x>12 & x<26
                %                fiber3=[ 1 0 0];
                %                  f3=1;
                % end








                %% geodesic C shape with fiber thickness one
                %
                % if( y>15 & y<32)
                % yy=(y-14);
                % xx=(x-16);
                %
                % %
                % if xx*xx+yy*yy>9*9 & xx*xx+yy*yy<10*10
                %     theta=atan(yy/xx);
                %     v=[sin(theta)  -cos(theta) 0];v=v/sqrt(v*v');
                %     fiber2=v;    f2=1;
                % end
                % end
                %% geodesic 3 Cubic not working

                %         xx=((16-x)^(-2));
                %         yy=abs(16-y);
                %         yy-xx^3;
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
%                 %%        Sine Shape  It works
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
 %%       parabola works
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
                %             m=0.4;%0.8 0.2 for one result sine and line 0.6 is merging
                %             if x-m*y>10 & x-m*y<15  % m=0.3,0.7,0.4--0.1,0.2 fails
                %                 theta=atan(m);
                %                 v=[ sin(theta) cos(theta) 0 ];%v=v/sqrt(v*v');
                %                 fiber2=v;
                %                 f2=1;
                %             end
                %             m=0.1; % pair (0.8,0.6) and (0.8 amd 0.7) 0.6 This one also
                %            if y-m*x>13 & y-m*x<19  % m=0.3,0.7,0.4--0.1,0.2 fails
                %            theta=atan(m);
                %            v=[ cos(theta) sin(theta) 0 ];
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
                %
                flagg=0;
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
                    pf3=0.0;
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
                    flagg=1;
                    fbckg_count=fbckg_count+1;
                    %                 pf1=0.8;%0.8,0.1,0.1 circle
                    %                 pf2=0.1;
                    %                 pf3=0.1;
                    %                fiber1=[1 0 0];
                    %                fiber2=[0 1 0];
                    %                fiber3=[0 0 1];

                end




                            if(flagg==0)
                S(x,y,1,i)=S(x,y,1,1)*(pf1*SimulateDWMRI(fiber1,GradientOrientations(i,:))+ ...
                    +pf2*SimulateDWMRI(fiber2,GradientOrientations(i,:))+...
                    +pf3*SimulateDWMRI(fiber3,GradientOrientations(i,:)));
                           else
                               S(x,y,1,i)=0.01;
                           end

                %             S(x,y,1,i)=S(x,y,1,1)*(pf1*SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+ pf2*SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:)));
                %            if ((x>13 && x<18) && (y>13 && y<18))
                %                S(x,y,1,i)=0.1;
                %            end
            end
        end
    end
    f1count=f1count/81;
    f2count=f2count/81;
    finter_count=finter_count/81;
    fbckg_count=fbckg_count/81;
    fregion=[f1count,f2count,finter_count,fbckg_count];
    %% Adding Riccian Noise
    % sigma=0.35;%0.02[0.01-0.09]in tha paper
    % for x=1:lsize
    %     for y=1:lsize
    %         for i=1:size(GradientOrientations,1)
    %             S(x,y,1,i)=sqrt( ( S(x,y,1,i)+sigma*randn(1) )^2+(sigma*randn(1))^2);
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


% A=squeeze(T6_CT(i,j,:,:));    
% tf = issymmetric(A);
% if(tf==1)
% d = eig(A)
% tol=length(d)*eps(max(d));
% isposdef = all(d > tol);
% issemidef = all(d > -tol);
% if(issemidef==1 || isposdef ==1 )
%  disp('Good SEMI SPD');   
% else
%     disp('Not SEMI SPD error');
%     f = msgbox('Not even semi SPD');
%     
% end
% end


        end
    end

    %% Till here We got the 6x6 tensors
    %%
    figure;
    delta=1;
    plotTensors(TensorCoefficients,delta,[321  1]);
    % figure;plotTensors(TensorCoefficients(:,16,16),1,[321 1] );
    figure;
    plotTensors(Tensor_ODF,delta,[321 1]);


    for i=1:size(S,1)

        for j=1:size(S,2)
            tensor1=squeeze(TensorCoefficients(:,i,j));
            W(i,j,:,:,:,:)=convertD2W(tensor1([15 5 1 12 3 10 11 8 7 14 13 9 4 6 2]));

            CT=squeeze(Tensor_ODF(:,i,j));
            W_CT(i,j,:,:,:,:)=convertD2W(CT([15 5 1 12 3 10 11 8 7 14 13 9 4 6 2]));

        end
    end
    for i=1:size(S,1)

        for j=1:size(S,2)

            T9fromFourthOrder(i,j,:,:)=[squeeze(W(i,j,:,:,1,1)), squeeze(W(i,j,:,:,1,2)), squeeze(W(i,j,:,:,1,3))
                squeeze(W(i,j,:,:,2,1)), squeeze(W(i,j,:,:,2,2)) , squeeze(W(i,j,:,:,2,3))
                squeeze(W(i,j,:,:,3,1)), squeeze(W(i,j,:,:,3,2)) , squeeze(W(i,j,:,:,3,3))];
            t11(:,:,i,j)=squeeze(W(i,j,:,:,1,1));
            t22(:,:,i,j)=squeeze(W(i,j,:,:,2,2));
            t33(:,:,i,j)=squeeze(W(i,j,:,:,3,3));
            %                                   t12(:,:,i,j)=squeeze(W(i,j,:,:,1,2))
% 
% A=squeeze(t33(:,:,i,j));    
% A=A-3;
% tf = issymmetric(A);
% if(tf==1)
% d = eig(A)
% tol=length(d)*eps(max(d));
% isposdef = all(d > tol);
% issemidef = all(d > -tol);
% if(issemidef==1 || isposdef ==1 )
%  disp('Good SEMI SPD');   
% else
%     disp('Not SEMI SPD error');
%     f = msgbox('Not even semi SPD');
%     
% end
% end


            
        end
    end
    %   t_sum=t1+t2+t3;
    for i=1:  size(T9D,2)
        for j=1:  size(T9D,1)
            %         disp('hello')

            T9_2D(i,j,:,:)= squeeze(t11(:,:,i,j)) + squeeze(t22(:,:,i,j));%+ squeeze(t33(:,:,i,j));
            %              T9_2D(i,j,:,:)=  squeeze(t11(:,:,i,j));
            %             T9_2D(i,j,:,:)= squeeze(t11(:,:,i,j))*squeeze(t22(:,:,i,j))*squeeze(t33(:,:,i,j));
            %
            %            if ((i>13 && i<18) && (j>13 && j<18))
            %                T9_2D(i,j,:,:)=T9_2D(i,j,:,:)/4;
            %            end

            Tensor2D(:,:,i,j)=T9_2D(i,j,:,:);
            chol( squeeze(T9_2D(i,j,:,:)) );
        end
    end
%      figure; plotDTI(Tensor2D,0.002);
    for i=1:  size(T9D,2)
        for j=1:  size(T9D,1)
            te=squeeze(T9D(i,j,:,:));
            [v, l]=eig(te);
            %            [v l]=svd(te);
            %         [ lambda, b]=max(diag(l));
            %         c=v(:,b);

            [ lambda, b]=sort(diag(l));

            c=v(:,(9));
            FstC=[c(1),c(2),c(3);c(4),c(5),c(6);c(7),c(8),c(9);].*lambda(9);
            c1=v(:,(8));
            ScndC=[c1(1),c1(2),c1(3);c1(4),c1(5),c1(6);c1(7),c1(8),c1(9);];%.*lambda(8);

            %          T9_2D(i,j,:,:)=FstC;
            % For display format Only
            %         c=v(:,9);
            t9_2D(:,:,i,j)=FstC;

            t_Eigen9_2D(i,j,:,:)=FstC;%[c(1),c(2),c(3);c(4),c(5),c(6);c(7),c(8),c(9);].*lambda;
            %           chol( squeeze( t_Eigen9_2D(i,j,:,:)) )
            %           T9_2D(i,j,:,:)=squeeze(t11(:,:,i,j))+ squeeze(t22(:,:,i,j))+ squeeze(t33(:,:,i,j))-t9_2D(:,:,i,j);
            chol( squeeze(T9_2D(i,j,:,:)) );
        end
    end
     end

    % Utility functions

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function W=convertD2W(D)
    %converts a vectorized version D([1:15]) with unique coeeficients
    %to a fully symmetric fourth order tensor W(:,:,:,:)
    global W;
    for i1=1:3
        for i2=1:3
            for i3=1:3
                for i4=1:3

                    ix=0;
                    iy=0;
                    iz=0;

                    if i1==1
                        ix=ix+1;
                    elseif i1==2
                        iy=iy+1;
                    elseif i1==3
                        iz=iz+1;
                    end

                    if i2==1
                        ix=ix+1;
                    elseif i2==2
                        iy=iy+1;
                    elseif i2==3
                        iz=iz+1;
                    end

                    if i3==1
                        ix=ix+1;
                    elseif i3==2
                        iy=iy+1;
                    elseif i3==3
                        iz=iz+1;
                    end

                    if i4==1
                        ix=ix+1;
                    elseif i4==2
                        iy=iy+1;
                    elseif i4==3
                        iz=iz+1;
                    end


                    if (ix==4)&(iy==0)&(iz==0)
                        W(i1,i2,i3,i4)=D(1)/computeFactor(4,0,0);
                    elseif (ix==0)&(iy==4)&(iz==0)
                        W(i1,i2,i3,i4)=D(2)/computeFactor(0,4,0);
                    elseif (ix==0)&(iy==0)&(iz==4)
                        W(i1,i2,i3,i4)=D(3)/computeFactor(0,0,4);
                    elseif (ix==2)&(iy==2)&(iz==0)
                        W(i1,i2,i3,i4)=D(4)/computeFactor(2,2,0);
                    elseif (ix==0)&(iy==2)&(iz==2)
                        W(i1,i2,i3,i4)=D(5)/computeFactor(0,2,2);
                    elseif (ix==2)&(iy==0)&(iz==2)
                        W(i1,i2,i3,i4)=D(6)/computeFactor(2,0,2);
                    elseif (ix==2)&(iy==1)&(iz==1)
                        W(i1,i2,i3,i4)=D(7)/computeFactor(2,1,1);
                    elseif (ix==1)&(iy==2)&(iz==1)
                        W(i1,i2,i3,i4)=D(8)/computeFactor(1,2,1);
                    elseif (ix==1)&(iy==1)&(iz==2)
                        W(i1,i2,i3,i4)=D(9)/computeFactor(1,1,2);
                    elseif (ix==3)&(iy==1)&(iz==0)
                        W(i1,i2,i3,i4)=D(10)/computeFactor(3,1,0);
                    elseif (ix==3)&(iy==0)&(iz==1)
                        W(i1,i2,i3,i4)=D(11)/computeFactor(3,0,1);
                    elseif (ix==1)&(iy==3)&(iz==0)
                        W(i1,i2,i3,i4)=D(12)/computeFactor(1,3,0);
                    elseif (ix==0)&(iy==3)&(iz==1)
                        W(i1,i2,i3,i4)=D(13)/computeFactor(0,3,1);
                    elseif (ix==1)&(iy==0)&(iz==3)
                        W(i1,i2,i3,i4)=D(14)/computeFactor(1,0,3);
                    elseif (ix==0)&(iy==1)&(iz==3)
                        W(i1,i2,i3,i4)=D(15)/computeFactor(0,1,3);
                    end
                end
            end
        end
    end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function counter=computeFactor(x,y,z)
    counter=0;
    for i1=1:3
        for i2=1:3
            for i3=1:3
                for i4=1:3
                    ix=0;
                    iy=0;
                    iz=0;
                    if i1==1
                        ix=ix+1;
                    elseif i1==2
                        iy=iy+1;
                    elseif i1==3
                        iz=iz+1;
                    end

                    if i2==1
                        ix=ix+1;
                    elseif i2==2
                        iy=iy+1;
                    elseif i2==3
                        iz=iz+1;
                    end

                    if i3==1
                        ix=ix+1;
                    elseif i3==2
                        iy=iy+1;
                    elseif i3==3
                        iz=iz+1;
                    end

                    if i4==1
                        ix=ix+1;
                    elseif i4==2
                        iy=iy+1;
                    elseif i4==3
                        iz=iz+1;
                    end
                    if (ix==x)&(iy==y)&(iz==z)
                        counter=counter+1;
                    end
                end
            end
        end
    end
    end
