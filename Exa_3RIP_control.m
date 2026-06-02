clear all; clc; close all;
syms x1 x2 x3 x4 u
format long;
%% control convexo del sistema pendulo de furuta de QUANSER
%Descriptor form:   E(x)*xp=A(x)*x+B(x)*u
%% Parámetros físicos (del paper)
g   = 9.81;

kt  = 0.00767;       % motor torque constant [N*m/A]
km  = 0.00767;       % back EMF constant [V/(rad/s)]
Kg  = 70;            % gear ratio
Rm  = 2.6;           % armature resistance [ohm]
eta_g = 0.9;         % gear efficiency
eta_m = 0.69;        % motor efficiency

Lp  = 0.3365;        % pendulum length [m]
Mp  = 0.127;         % pendulum mass [kg]
Jp  = 0.0012;        % pendulum inertia [kg*m^2]
Bp  = 0.0024;        % pendulum damping

Lr  = 0.2159;        % rotary arm length [m]
Mr  = 0.257;         % rotary arm mass [kg]
Jr  = 9.9829e-4;     % rotary arm inertia [kg*m^2]
Br  = 0.0024;        % rotary arm damping

T1 = Rm/(eta_g*Kg*eta_m*kt);
T2 = Mp*Lp^2;
T3 = Mp*Lr^2;
T4 = Mp*Lp*Lr;
T5 = Kg*km+T1*Br;
T6 = Jp+0.25*T2;

%%% decay rate
alfa=0.4; %%control
%% Modelado convexo
% Vector de premisas: no linealidades
%medibles
z(1) = sin(x2);
z(2) = cos(x2);
z(3) = sin(x2)/x2;

%No medibles
tzeta(1) = x3;
tzeta(2) = x4;
tzeta(3)=sin(x2)^2;

%region de modelado
x1_limit=deg2rad(90); %limit on the arm
x2_limit=deg2rad(14); %limit  for the pendulum
x3_limit=5; %limit of arm's velocity
x4_limit=10; %limit of pendulum's velocity

x1=-x1_limit:0.0001:x1_limit;
x2=-x2_limit:0.0001:x2_limit;
x3=-x3_limit:0.0001:x3_limit;
x4=-x4_limit:0.0001:x4_limit;


%%% cotas
for i=1:numel(z)
    z_bound(i,1)=min(eval(z(i)));
    z_bound(i,2)=max(eval(z(i)));
end

for i=1:numel(tzeta)
    tzeta_bound(i,1)=min(eval(tzeta(i)));
    tzeta_bound(i,2)=max(eval(tzeta(i)));
end

% Pesos
for i = 1:numel(z)
    w(i,1) = (z_bound(i,2)-z(i)) / (z_bound(i,2) - z_bound(i,1));
    w(i,2) = 1 - w(i,1);
end


% Funciones de membresía (5 premisas => 2^5 = 32)
h(1) = w(1,1)*w(2,1)*w(3,1);
h(2) = w(1,1)*w(2,1)*w(3,2);
h(3) = w(1,1)*w(2,2)*w(3,1);
h(4) = w(1,1)*w(2,2)*w(3,2);

h(5) = w(1,2)*w(2,1)*w(3,1);
h(6) = w(1,2)*w(2,1)*w(3,2);
h(7) = w(1,2)*w(2,2)*w(3,1);
h(8) = w(1,2)*w(2,2)*w(3,2);


r = 2^numel(z); %cantidad de vertices a generar
rho=2^numel(tzeta);
%%% Matrices A  B y E
A = cell(r,1); B = cell(r,1); E = cell(r,1);
i=0;
for i1=1:2
    for i2=1:2
        for i3=1:2
            
            i = i + 1;
            j=0;
            for j1=1:2
                for j2=1:2
                    for j3=1:2
                        j=j+1;
                        
                        E{i,j} = [1, 0, 0,0;...
                            0, 1, 0, 0;...
                            0, 0, T1*(0.25*T2*(1-tzeta_bound(3,j3))+Jr+T3), -0.5*T1*T4*z_bound(2,i2);...
                            0, 0, -0.5*T4*z_bound(2,i2), T6];
                        
                        A{i,j} =[0, 0, 1, 0;...
                            0, 0, 0, 1;...
                            0, 0, -0.5*T1*T2*z_bound(1,i1)*z_bound(2,i2)*tzeta_bound(2,j2)-T5, -0.5*T1*T4*z_bound(1,i1)*tzeta_bound(2,j2);...
                            0, 0.5*Mp*Lp*g*z_bound(3,i3), 0.25*T2*z_bound(1,i1)*z_bound(2,i2)*tzeta_bound(1,j1), -Bp];
                        
                        B = [0;0;1;0];
                        
                    end
                end
            end
        end
    end
end
%% Calculo de ley de control

[n, m] = size(B);
setlmis([]); %%%inicializa un grupo de LMIs

%decision variables
P1=lmivar(1,[n 1]);
P3=lmivar(2,[n n]);
P4=lmivar(2,[n n]);
for j=1:r
    M{j}=lmivar(2,[m n]);
end

%%% LMIs for control
LMI_ctr1=newlmi;
lmiterm([-LMI_ctr1 1 1 P1],1,1);
for i=1:r
    for j=1:rho
        LMI_ctr2=newlmi;
        lmiterm([LMI_ctr2 1 1 P3],1,1,'s');
        lmiterm([LMI_ctr2 1 1 P1],1,2*alfa);
        lmiterm([LMI_ctr2 2 1 P1],A{i,j},1);
        lmiterm([LMI_ctr2 2 1 M{i}],B,1);
        lmiterm([LMI_ctr2 2 1 P3],-E{i,j},1);
        lmiterm([LMI_ctr2 2 1 -P4],1,1);
        lmiterm([LMI_ctr2 2 2 P4],-E{i,j},1,'s');
    end
end
%%% LMIs restriccion a la entrada
x0 = [deg2rad(45);deg2rad(15);1;3]; %condiciones iniciales en las que se pondra en el pendulo
mu = 10; %voltaje de restriccion

%%% LMI restriccion a la entrada (1)
LMI_input_1 = newlmi;
% (1,1) = 1
lmiterm([-LMI_input_1 1 1 0], 1);
% (2,1) = x0
lmiterm([-LMI_input_1 2 1 0], x0);
% (2,2) = P1
lmiterm([-LMI_input_1 2 2 P1], 1, 1);


%%% LMI restriccion a la entrada (2)
for i = 1:r
    LMI_input_2 = newlmi;
    % (1,1) = P1
    lmiterm([-LMI_input_2 1 1 P1], 1, 1);
    % (2,1) = M_i
    lmiterm([-LMI_input_2 2 1 M{i}], 1, 1);
    % (2,2) = mu^2 * I_m
    lmiterm([-LMI_input_2 2 2 0], mu^2);
end


ejemplo=getlmis;
[tmin, xfeas]=feasp(ejemplo); %%manda a llamar el solver
K_block = zeros(r,4);
if tmin<0
    disp('FACTIBLE');
    P1=dec2mat(ejemplo,xfeas,P1);
    P3=dec2mat(ejemplo,xfeas,P3);
    P4=dec2mat(ejemplo,xfeas,P4);
    for j=1:r
        M{j}=dec2mat(ejemplo,xfeas,M{j});
        K{j}=M{j}*inv(P1);
        K_block(j,:) = K{j};
    end
    fprintf('K block:\n');
    disp(K_block);
else
    disp('No FACTIBLE');
end
