function[All_Tensor_Coff]=Fibers_0_to_5_specific_iteration()
S0=1.0;
UnitVectors ;
n_grad=81;
GradientOrientations=[1 0 0; g(1:n_grad,:)];
b_value=[10;ones(n_grad,1)*1500];
bvalue=1500;
f=sqrt(2);
iphi=0;
ith=0;
th=90;
phi=0;
th_inc=0;
phi_inc=0;
diag_sz=.33;
inc_diag_sz=0;
% Adding Riccian Noise
%    sigma=0.09;%0.02[0.01-0.09]in tha paper

f1=[1,0,0]; f2=[0,1,0]; f3=[0,0,1];
max_fiber_per_voxel=3;
no_of_itration=6;
All_Tensor_Coff=zeros(15,max_fiber_per_voxel+1,no_of_itration);
for ii=1 : no_of_itration
    phi=phi+phi_inc;
    %   th=th+th_inc;
    
    
    fiber_orientation0=[0.0001 0.0001  0.0001];
    
    fiber_orientation1=[cos(phi*pi/180)*sin(th*pi/180) sin(phi*pi/180)*sin(th*pi/180) cos(th*pi/180)];
    % th=10*randi([-9,9]);  phi=10*randi([-18,18]);
    % fiber_orientation1=[.99 .09985 .09985];
    fiber_orientation2=[cos((phi+90)*pi/180)*sin(th*pi/180) sin((phi+90)*pi/180)*sin(th*pi/180) cos(th*pi/180)];
    % th=10*randi([-9,9]);  phi=10*randi([-18,18]);
    % fiber_orientation2=[.09985 .99 .09985];
    %  fiber_orientation3=[cos((phi+90)*pi/180)*sin((th+90)*pi/180) sin((phi+90)*pi/180)*sin((th+90)*pi/180) cos((th+90)*pi/180)];
    fiber_orientation3=[0 0 1];
    % th=10*randi([-9,9]);  phi=10*randi([-18,18]);
    % fiber_orientation3=[.09985 .09985 .99];
    % fiber_orientation4=[cos(phi*pi/180)*sin(th*pi/180) sin(phi*pi/180)*sin(th*pi/180) cos(th*pi/180)];
    fiber_0=zeros(size(GradientOrientations,1),1);
    fiber_1 = zeros(size(GradientOrientations,1),1);
    fiber_1_and_2 = zeros(size(GradientOrientations,1),1);
    fiber_1_2_and_3 = zeros(size(GradientOrientations,1),1);
    % fiber_1_2_3_and_4 = zeros(size(GradientOrientations,1),1);
    % std_fiber_xyz=zeros(size(GradientOrientations,1),1);
    % diag_sz=diag_sz+inc_diag_sz;
    
    for i=1:size(GradientOrientations,1)
        % Fiber 0 / no fiber at voxel only symmetric/isotropic tensor.
        fiber_0(i)=diag_sz;
        S(1,1,i)=diag_sz;
        % Fiber 1 / only single fiber at one voxel.
        fiber_1(i)=S0* (SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:)));
        S(1,2,i)=S0* (SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:)));
        % Fiber 1 and 2 / Two crossing fibers at one voxel.
        fiber_1_and_2(i)=S0* (SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+ SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:)))/2;
        S(1,3,i)=S0* (SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+ SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:)))/2;
        % Fiber 1,2 and 3 / Three crossing fibers at one voxel.
        fiber_1_2_and_3(i)=S0* (SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+ SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:))+SimulateDWMRI(fiber_orientation3,GradientOrientations(i,:)))/3;
        S(1,4,i)=S0* (SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+ SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:))+SimulateDWMRI(fiber_orientation3,GradientOrientations(i,:)))/3;
        
        %
        %    % Fiber 1,2,3 and 4 / Three crossing fibers at one voxel.
        %    fiber_1_2_3_and_4(i)=S0* (SimulateDWMRI(fiber_orientation1,GradientOrientations(i,:))+ SimulateDWMRI(fiber_orientation2,GradientOrientations(i,:))+SimulateDWMRI(fiber_orientation3,GradientOrientations(i,:))+SimulateDWMRI(fiber_orientation4,GradientOrientations(i,:)))/4;
        %
        %    % Std fiber in x y z direction / Three crossing fibers at one voxel at 90 each.
        %    std_fiber_xyz(i)=S0* (SimulateDWMRI(f1,GradientOrientations(i,:))+ SimulateDWMRI(f2,GradientOrientations(i,:))+SimulateDWMRI(f3,GradientOrientations(i,:)))/3;
    end
    % for i=1:size(GradientOrientations,1)
    % %% add richian noise
    %  fiber_0(i)=sqrt( ( fiber_0(i)+sigma*randn(1) )^2+(sigma*randn(1))^2);
    % fiber_1(i)=sqrt( ( fiber_1(i)+sigma*randn(1) )^2+(sigma*randn(1))^2);
    % fiber_1_and_2(i)=sqrt( ( fiber_1_and_2(i)+sigma*randn(1) )^2+(sigma*randn(1))^2);
    % fiber_1_2_and_3(i)=sqrt( ( fiber_1_2_and_3(i)+sigma*randn(1) )^2+(sigma*randn(1))^2);
    % % S(x,y,1,i)=sqrt((S(x,y,1,i)+sigma*randn(1))^2+(sigma*randn(1))^2);
    %
    % end
    
    
    
    order=4;
    G=constructMatrixOfMonomials(GradientOrientations, order);
    C=constructSetOf321Polynomials(order)';
    P=G*C;
    P=[-diag(b_value)*P ones(size(GradientOrientations,1),1)];
    BG=constructMatrixOfIntegrals(GradientOrientations, order, 100);
    B=BG*C;
    y_fiber_0=squeeze(log(fiber_0));
    x_fiber_0=lsqnonneg(P, y_fiber_0);
    y_fiber_1=squeeze(log(fiber_1));
    x_fiber_1=lsqnonneg(P, y_fiber_1);
    y_fiber_1_and_2=squeeze(log(fiber_1_and_2));
    x_fiber_1_and_2= lsqnonneg(P, y_fiber_1_and_2);
    y_fiber_1_2_and_3=squeeze(log(fiber_1_2_and_3));
    x_fiber_1_2_and_3= lsqnonneg(P, y_fiber_1_2_and_3);
    % x_fiber_1_2_3_and_4= lsqnonneg(P, fiber_1_2_3_and_4/S0);
    % x_std_fiber_xyz= lsqnonneg(P, std_fiber_xyz/S0);
    
    Old_Tensor_Coeff_fiber_0(:,ii) = C * x_fiber_0([1:321]);
    Old_Tensor_Coeff_fiber_1(:,ii) = C * x_fiber_1([1:321]);
    Old_Tensor_Coeff_fiber_1_and_2(:,ii) = C * x_fiber_1_and_2([1:321]);
    Old_Tensor_Coeff_fiber_1_2_and_3(:,ii) = C * x_fiber_1_2_and_3([1:321]);
    %via CT-fod
    x1_fiber_0=lsqnonneg(B, exp(-bvalue*G*Old_Tensor_Coeff_fiber_0(:,ii)));
    New_Tensor_Coeff_fiber_0(:,ii) = C *  x1_fiber_0;
    
    x1_fiber_1=lsqnonneg(B, exp(-bvalue*G*Old_Tensor_Coeff_fiber_1(:,ii)));
    New_Tensor_Coeff_fiber_1(:,ii) = C *  x1_fiber_1;
    
    x1_fiber_1_and_2=lsqnonneg(B, exp(-bvalue*G*Old_Tensor_Coeff_fiber_1_and_2(:,ii)));
    New_Tensor_Coeff_fiber_1_and_2(:,ii) = C *  x1_fiber_1_and_2;
    
    x1_fiber_1_2_and_3=lsqnonneg(B, exp(-bvalue*G*Old_Tensor_Coeff_fiber_1_2_and_3(:,ii)));
    New_Tensor_Coeff_fiber_1_2_and_3(:,ii) = C *  x1_fiber_1_2_and_3;
    
    OldD=zeros(15,4);
    OldD(:,1)=Old_Tensor_Coeff_fiber_0(:,ii);
    OldD(:,2)=Old_Tensor_Coeff_fiber_1(:,ii);
    OldD(:,3)=Old_Tensor_Coeff_fiber_1_and_2(:,ii);
    OldD(:,4)=Old_Tensor_Coeff_fiber_1_2_and_3(:,ii);
    
    D=zeros(15,4);
    D(:,1)=New_Tensor_Coeff_fiber_0(:,ii);
    D(:,2)=New_Tensor_Coeff_fiber_1(:,ii);
    D(:,3)=New_Tensor_Coeff_fiber_1_and_2(:,ii);
    D(:,4)=New_Tensor_Coeff_fiber_1_2_and_3(:,ii);
    % Old_Tensor_Coeff_fiber_1_2_3_and_4 = C * x_fiber_1_2_3_and_4;
    % Old_Tensor_Coeff_std_fiber_xyz = C * x_std_fiber_xyz;
    % figure; plotTensors(Old_Tensor_Coeff_fiber_0, .002);
    % figure; plotTensors(Old_Tensor_Coeff_fiber_1, .002);
    % figure; plotTensors(Old_Tensor_Coeff_fiber_1_and_2, .002);
    % figure; plotTensors(Old_Tensor_Coeff_fiber_1_2_and_3, .002);
    % % figure; plotTensors(Old_Tensor_Coeff_fiber_1_2_3_and_4, .002);
    % % figure; plotTensors(Old_Tensor_Coeff_std_fiber_xyz, .002);
    
    All_Tensor_Coff(:,1,ii)=D(:,2);
    All_Tensor_Coff(:,2,ii)=D(:,2);
    All_Tensor_Coff(:,3,ii)=D(:,2);
    All_Tensor_Coff(:,4,ii)=D(:,3);
    All_Tensor_Coff(:,5,ii)=D(:,3);
    All_Tensor_Coff(:,6,ii)=D(:,3);
    All_Tensor_Coff(:,7,ii)=D(:,4);
    All_Tensor_Coff(:,8,ii)=D(:,4);
    All_Tensor_Coff(:,9,ii)=D(:,4);
    All_Tensor_Coff(:,10,ii)=D(:,3);
    All_Tensor_Coff(:,11,ii)=D(:,3);
    All_Tensor_Coff(:,12,ii)=D(:,3);
    All_Tensor_Coff(:,13,ii)=D(:,2);
    All_Tensor_Coff(:,14,ii)=D(:,2);
    All_Tensor_Coff(:,15,ii)=D(:,2);
    All_Tensor_Coff(:,16,ii)=D(:,4);
    All_Tensor_Coff(:,17,ii)=D(:,4);
    All_Tensor_Coff(:,18,ii)=D(:,4);
    All_Tensor_Coff(:,19,ii)=D(:,2);
    All_Tensor_Coff(:,20,ii)=D(:,2);
    All_Tensor_Coff(:,21,ii)=D(:,2);  
    
    Old_All_Tensor_Coff(:,1,ii)=OldD(:,1);
    Old_All_Tensor_Coff(:,2,ii)=OldD(:,2);
    Old_All_Tensor_Coff(:,3,ii)=OldD(:,3);
    Old_All_Tensor_Coff(:,4,ii)=OldD(:,4);
end

% figure;
% plotTensors(Old_All_Tensor_Coff,1,[321  1]);

figure;
plotTensors(All_Tensor_Coff,1,[321  1]);



end