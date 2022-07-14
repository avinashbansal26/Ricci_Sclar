 function[Tensor_ODF]=Fibers_specific_iteration()
S0=1.0;
UnitVectors ;
n_grad=81;
GradientOrientations=[1 0 0; g(1:n_grad,:)];
b_value=[10;ones(n_grad,1)*1500];
bvalue=1500;
f1=[1,0,0]; f2=[0,1,0]; f3=[0,0,1];
max_fiber_per_voxel=3;
no_of_itration=4;
% S=ones(no_of_itration,max_fiber_per_voxel+1,size(GradientOrientations,1));
All_Tensor_Coff=zeros(15,max_fiber_per_voxel+1,no_of_itration);   
for ii=1 : no_of_itration
    thx=rand*pi/3;
    thy=rand*pi/4;
    thz=rand*pi/2;
    rot_x=[1 0 0; 0 cos(thx) -sin(thx); 0 sin(thx) cos(thx)];
    rot_y=[cos(thy) 0 sin(thy); 0 1 0; -sin(thy) 0 cos(thy)];
    rot_z=[cos(thz) -sin(thz) 0; sin(thz) cos(thz) 0; 0 0 1];
    R=rot_z*rot_y*rot_x;
    fiber1=R*f1';
    fiber2=R*f2';
    fiber3=R*f3';
    for i=1:size(GradientOrientations,1)
        S(ii,1,i)=.001;
        S(ii,2,i)=S0* (SimulateDWMRI(fiber1,GradientOrientations(i,:)));
        S(ii,3,i)=S0* (SimulateDWMRI(fiber1,GradientOrientations(i,:))+ SimulateDWMRI(fiber2,GradientOrientations(i,:)))/2;
        S(ii,4,i)=S0* (SimulateDWMRI(fiber1,GradientOrientations(i,:))+ SimulateDWMRI(fiber2,GradientOrientations(i,:))+SimulateDWMRI(fiber3,GradientOrientations(i,:)))/3;
    end
end

    order=4;
    G=constructMatrixOfMonomials(GradientOrientations, order);
    C=constructSetOf321Polynomials(order)';
    P=G*C;
    P=[-diag(b_value)*P ones(size(GradientOrientations,1),1)];
    BG=constructMatrixOfIntegrals(GradientOrientations, order, 100);
    B=BG*C;
    for i=1:size(S,1)
        for j=1:size(S,2)
            y=squeeze(log(S(i,j,:)));
            x=lsqnonneg(P, y);
            UniqueTensorCoefficients = C * x(1:321);
            TensorCoefficients(:,i,j)=UniqueTensorCoefficients;
            D=UniqueTensorCoefficients;
            %% For CT ODF
            x1=lsqnonneg(B, exp(-bvalue*G*D));
            CT_Coeff = C * x1;
            Tensor_ODF(:,i,j)=CT_Coeff;            
        end
    end

figure;
plotTensors(TensorCoefficients,1,[321  1]);
figure;
plotTensors(Tensor_ODF,1,[321  1]);
end