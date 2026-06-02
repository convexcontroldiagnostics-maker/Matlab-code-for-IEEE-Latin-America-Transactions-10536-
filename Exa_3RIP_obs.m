clear all; clc; close all;
syms x1 x2 x3 x4 x3g x4g u
%%% Observador convexo del sistema pendulo de furuta de QUANSER
%% Parámetros físicos (del paper Lendek)
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
beta=0.1; %%Observer

%% Modelado convexo
% Vector de premisas
%Conocidas
tz(1) = sin(x2);
tz(2) = cos(x2);
tz(3) = x3g;
tz(4) = x4g;
tz(5) = (cos(x2))^2;
%No conocidas
tzeta(1) = x3;
tzeta(2) = x4;

%region de modelado
x1_limit=deg2rad(90);
x2_limit=deg2rad(15);
x3_limit=4;
x4_limit=4;

x1=-x1_limit:0.0001:x1_limit;
x2=-x2_limit:0.0001:x2_limit;
x3=-x3_limit:0.0001:x3_limit;
x4=-x4_limit:0.0001:x4_limit;
x3g=-x3_limit:0.0001:x3_limit;
x4g=-x4_limit:0.0001:x4_limit;


%%% Modeling for the observer
for i=1:numel(tz)
    tz_bound(i,1)=min(eval(tz(i)));
    tz_bound(i,2)=max(eval(tz(i)));
end

for i=1:numel(tzeta)
    tzeta_bound(i,1)=min(eval(tzeta(i)));
    tzeta_bound(i,2)=max(eval(tzeta(i)));
end

re=2^(numel(tz)); %%% for the observer
rho=2^(numel(tzeta)); %%% for the observer

%%% vertex matrices:
j=0;
for j1=1:2
    for j2=1:2
        j=j+1;i=0;
        for i1=1:2
            for i2=1:2
                for i3=1:2
                    for i4=1:2
                        for i5=1:2
                            
                            i=i+1;
                            Ee{i,j} = [1, 0, 0, 0;...
                                0, 1, 0, 0;...
                                0, 0, T1*(0.25*T2*(1-tz_bound(2,i2)*tz_bound(5,i5))+Jr+T3), -0.5*T1*T4*tz_bound(2,i2);
                                0, 0, -0.5*T4*tz_bound(2,i2), T6];
                            
                            
                            Ae{i,j}= [0, 0, 1, 0;
                                0, 0, 0, 1;
                                0, 0, -0.5*T1*T2*tz_bound(1,i1)*tz_bound(2,i2)*tz_bound(4,i4), -0.5*T1*T2*tz_bound(1,i1)*(tz_bound(2,i2)*tzeta_bound(1,i1)+tzeta_bound(2,i2)+tz_bound(4,i4));
                                0, 0, 0.25*T2*tz_bound(2,i2)*tz_bound(1,i1)*(tzeta_bound(1,i1)+tz_bound(3,i3)), -Bp];
                            
                        end
                    end
                end
            end
        end
    end
end
Ce = [1, 0, 0, 0;0, 1, 0, 0];
%% Calculo de ley de control para observador

[q, n]=size(Ce);
setlmis([]); %%%inicializa un grupo de LMIs

P1e=lmivar(1,[n 1]);
P3e=lmivar(2,[n n]);
P4e=lmivar(2,[n n]);
for j=1:re
    M1e{j}=lmivar(2,[n q]);
    M2e{j}=lmivar(2,[n q]);
end

LMI_obs1=newlmi;
lmiterm([-LMI_obs1 1 1 P1e],1,1);

for i=1:re
    for k=1:rho
        LMI_obs2=newlmi;
        lmiterm([LMI_obs2 1 1 -P3e],1,Ae{i,k},'s');
        lmiterm([LMI_obs2 1 1 M1e{i}],1,-Ce,'s');
        lmiterm([LMI_obs2 1 1 P1e],1,2*beta);
        lmiterm([LMI_obs2 2 1 -P4e],1,Ae{i,k});
        lmiterm([LMI_obs2 2 1 M2e{i}],1,-Ce);
        lmiterm([LMI_obs2 2 1 P1e],1,1);
        lmiterm([LMI_obs2 2 1 P3e],-Ee{i}',1);
        lmiterm([LMI_obs2 2 2 -P4e],1,-Ee{i},'s');
    end
end

mu=0.02;
H1=lmivar(1,[n 1]);
H2=lmivar(1,[n 1]);
for i=1:re
    LMI_bound=newlmi;
    lmiterm([-LMI_bound 1 1 0],mu);
    lmiterm([-LMI_bound 2 1 M1e{i}],1,1);
    lmiterm([-LMI_bound 3 1 M2e{i}],1,1);
    lmiterm([-LMI_bound 2 2 P1e],2,1);
    lmiterm([-LMI_bound 2 2 H1],-1,1);
    lmiterm([-LMI_bound 3 2 P3e],1,1);
    lmiterm([-LMI_bound 3 3 P4e],1,1,'s');
    lmiterm([-LMI_bound 3 3 H2],-1,1);
end

ejemplo=getlmis;
[tmin, xfeas]=feasp(ejemplo); %%manda a llamar el solver

M1e_vector = zeros(4,2);
M2e_vector = zeros(4,2);

if tmin<0
    disp('FACTIBLE');
    format long;
    P1e=dec2mat(ejemplo,xfeas,P1e);
    P3e=dec2mat(ejemplo,xfeas,P3e);
    P4e=dec2mat(ejemplo,xfeas,P4e);
    
    fprintf('Pbar = inv(transpose([P1e, 0;P3e, P4e])):\n');
    Pbar = inv(transpose([P1e, zeros(4,4);P3e, P4e]));
    disp(Pbar);
    
    
    for j=1:re
        
        M1e{j}=dec2mat(ejemplo,xfeas,M1e{j});
        M2e{j}=dec2mat(ejemplo,xfeas,M2e{j});
        
        cols = (1:2) + (j-1)*2;
        M1e_vector(:, cols) = M1e{j};
        M2e_vector(:, cols) = M2e{j};
    end
    
    fprintf('M1e_vector:\n');
    disp(M1e_vector);
    fprintf('M2e_vector:\n');
    disp(M2e_vector);
    
    
else
    disp('No FACTIBLE');
end
