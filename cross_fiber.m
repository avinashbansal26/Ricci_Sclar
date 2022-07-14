 function [Tensor_ODF]=cross_fiber()
order=4;%In this demo we compute 4th-order tensorsfor ii = 1:9
S0=1.0;
UnitVectors;
n_grad=81;
GradientOrientations=[1 0 0;g([1:n_grad],:)];
b_value=[10;ones(n_grad,1)*1500];
lsize=16;
bvalue=1500;
S=ones(lsize,lsize,1,size(GradientOrientations,1));
fiber1=[1 0 0]; fiber2=[0 1 0];
for x=1:lsize
    for y=1:lsize      
        for i=2:size(GradientOrientations,1)
            %%  cross fiber Code you can use this code   
            if y<=11 && y>=6
               f1=1; 
            else
                f1=0;
            end          
           %Define the underlying fiber orientation of fiber bundle #2 (i.e. we have 2-fiber crossing)
           if x<=11 && x>=6
               f2=1;
           else
               f2=0;
           end
                S(x,y,i)=(f1*SimulateDWMRI(fiber1,GradientOrientations(i,:))+ ...
                    f2*SimulateDWMRI(fiber2,GradientOrientations(i,:)))/2;
                if  S(x,y,1,i)==0
                     S(x,y,1,i)=.01;
                end
        end
%         %% Adding Riccian Noise
%         sigma=0.05;%0.02[0.01-0.09]in tha paper
%         for i=1:size(GradientOrientations,1)
%             S(x,y,1,i)=sqrt((S(x,y,1,i)+sigma*randn(1))^2+(sigma*randn(1))^2);    
%         end
    end
end
   
%Construct all possible monomials of a specific order
G=constructMatrixOfMonomials(GradientOrientations, order); %computes G from section 5.1 (ISBI'10)
%Construct set of polynomial coefficients C
C=constructSetOf321Polynomials(order)'; %computes C from section 5.1 (ISBI'10)
P=G*C;
P=[-diag(b_value)*P ones(size(GradientOrientations,1),1)];
%%For CT ODF
BG=constructMatrixOfIntegrals(GradientOrientations, order, 100);
B=BG*C;
for i=1:size(S,1)
    for j=1:size(S,2)
        y=squeeze(log(S(i,j,:)));
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
    end
end
delta=1;
scaling=1;
visibilty=-1;
figure; plotTensors(TensorCoefficients,delta,[321 scaling visibilty]);
figure; plotTensors(Tensor_ODF,delta,[321 scaling visibilty]);
 end
