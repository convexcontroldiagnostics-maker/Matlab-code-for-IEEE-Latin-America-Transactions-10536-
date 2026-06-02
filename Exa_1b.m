%%% LMI conditions in Theorem 1 and 2 for observer design
%%% Example 1

clear all; clc;


global P1 P3 P4 K z_bound P1e P3e P4e M1e M2e tz_bound tzeta_bound

syms x1 x2 x3 x4 x1g x2g x3g x4g

%%% decay rate
alfa=0.5; %%control
beta=1; %%Observer
mu=5; %% input bound
x0=[-1.2 0.5 0.7 -0.6]';
%%% Modeling region:
x1_limit=1.4;
x3_limit=0.7;

%%% Premise vector
z(1)=x1^2;
z(2)=sin(x3)/x3;

tz(1)=x1^2;


tzeta(1)=x3^2+x3*x3g+x3g^2;




x1=-x1_limit:0.0001:x1_limit;
x3=-x3_limit:0.0001:x3_limit;
x3g=-x3_limit:0.0001:x3_limit;


%%% Modeling for the controller


for i=1:numel(z)
    z_bound(i,1)=min(eval(z(i)));
    z_bound(i,2)=max(eval(z(i)));
end
%%% vertex matrices:
i=0;
for i1=1:2
    for i2=1:2
        i=i+1;
        E{i}=eye(4);
        A{i}=[1 1 z_bound(2,i2) -0.1;
            1 -2 0 0;
            1 z_bound(1,i1) -0.3 0;
            0 0 z_bound(2,i2) -1];
        B{i}=[z_bound(1,i1)+1;0;0;0];
    end
end

[n m]=size(B{1});



%%% Modeling for the observer


for i=1:numel(tz)
    tz_bound(i,1)=min(eval(tz(i)));
    tz_bound(i,2)=max(eval(tz(i)));
end

for i=1:numel(tzeta)
    tzeta_bound(i,1)=min(eval(tzeta(i)));
    tzeta_bound(i,2)=max(eval(tzeta(i)));
end

%%% vertex matrices:
j=0;
for j1=1:2
    j=j+1;i=0;
    for i1=1:2
        i=i+1;
        Ee{i}=eye(4);
        Ae{i,j}=[1 1 1-1/6*tzeta_bound(1,j1) -0.1;
            1 -2 0 0;
            1 tz_bound(1,i1) -0.3 0;
            0 0 1-1/6*tzeta_bound(1,j1) -1];
        Ce{i,j}=[0 1 0 tz_bound(1,i1)+1;
            1 0 0 0];
        
    end
end


[q n]=size(Ce{1,1});

factible=0;

r=2^(numel(z)); %%% for the controller
re=2^(numel(tz)); %%% for the observer
rho=2^(numel(tzeta)); %%% for the observer


setlmis([]); %%%inicializa un grupo de LMIs

%decision variables
P1=lmivar(1,[n 1]);
for j=1:r
    P3{j}=lmivar(2,[n n]);
    P4{j}=lmivar(2,[n n]);
    M{j}=lmivar(2,[m n]);
end

P1e=lmivar(1,[n 1]);

for j=1:re
    M1e{j}=lmivar(2,[n q]);
    M2e{j}=lmivar(2,[n q]);
    P3e{j}=lmivar(2,[n n]);
P4e{j}=lmivar(2,[n n]);
end


%%% LMIs for control
LMI_ctr1=newlmi;
lmiterm([-LMI_ctr1 1 1 P1],1,1);
for i=1:r
    for j=1:r
        LMI_ctr2=newlmi;
        lmiterm([LMI_ctr2 1 1 P3{i}],1,2/(r-1),'s');
        lmiterm([LMI_ctr2 1 1 P1],2/(r-1),2*alfa);
        lmiterm([LMI_ctr2 2 1 P1],A{i},2/(r-1));
        lmiterm([LMI_ctr2 2 1 M{i}],B{i},2/(r-1));
        lmiterm([LMI_ctr2 2 1 P3{i}],-E{i},2/(r-1));
        lmiterm([LMI_ctr2 2 1 -P4{i}],1,2/(r-1));
        lmiterm([LMI_ctr2 2 2 P4{i}],-E{i},2/(r-1),'s');
        
        lmiterm([LMI_ctr2 1 1 P3{j}],1,1,'s');
        lmiterm([LMI_ctr2 1 1 P1],1,2*alfa);
        lmiterm([LMI_ctr2 2 1 P1],A{i},1);
        lmiterm([LMI_ctr2 2 1 M{j}],B{i},1);
        lmiterm([LMI_ctr2 2 1 P3{j}],-E{i},1);
        lmiterm([LMI_ctr2 2 1 -P4{j}],1,1);
        lmiterm([LMI_ctr2 2 2 P4{j}],-E{i},1,'s');
        
        lmiterm([LMI_ctr2 1 1 P3{i}],1,1,'s');
        lmiterm([LMI_ctr2 1 1 P1],1,2*alfa);
        lmiterm([LMI_ctr2 2 1 P1],A{j},1);
        lmiterm([LMI_ctr2 2 1 M{i}],B{j},1);
        lmiterm([LMI_ctr2 2 1 P3{i}],-E{j},1);
        lmiterm([LMI_ctr2 2 1 -P4{i}],1,1);
        lmiterm([LMI_ctr2 2 2 P4{i}],-E{j},1,'s');
        
    end
end


LMI0=newlmi;
lmiterm([-LMI0 1 1 0],1);
lmiterm([-LMI0 2 1 0],x0);
lmiterm([-LMI0 2 2 P1],1,1);

for j=1:r
    LMI0=newlmi;
    lmiterm([-LMI0 1 1 P1],1,1);
    lmiterm([-LMI0 2 1 M{j}],1,1);
    lmiterm([-LMI0 2 2 0],mu^2);
end

LMI_obs1=newlmi;
lmiterm([-LMI_obs1 1 1 P1e],1,1);

for i=1:re
    for j=1:re
        for k=1:rho
            LMI_obs2=newlmi;
            lmiterm([LMI_obs2 1 1 -P3e{i}],2/(r-1),Ae{i,k},'s');
            lmiterm([LMI_obs2 1 1 M1e{i}],2/(r-1),-Ce{i,k},'s');
            lmiterm([LMI_obs2 1 1 P1e],2/(r-1),2*beta);
            lmiterm([LMI_obs2 2 1 -P4e{i}],2/(r-1),Ae{i,k});
            lmiterm([LMI_obs2 2 1 M2e{i}],2/(r-1),-Ce{i,k});
            lmiterm([LMI_obs2 2 1 P1e],2/(r-1),1);
            lmiterm([LMI_obs2 2 1 P3e{i}],-Ee{i}',2/(r-1));
            lmiterm([LMI_obs2 2 2 -P4e{i}],2/(r-1),-Ee{i},'s');
            
            lmiterm([LMI_obs2 1 1 -P3e{j}],1,Ae{i,k},'s');
            lmiterm([LMI_obs2 1 1 M1e{j}],1,-Ce{i,k},'s');
            lmiterm([LMI_obs2 1 1 P1e],1,2*beta);
            lmiterm([LMI_obs2 2 1 -P4e{j}],1,Ae{i,k});
            lmiterm([LMI_obs2 2 1 M2e{j}],1,-Ce{i,k});
            lmiterm([LMI_obs2 2 1 P1e],1,1);
            lmiterm([LMI_obs2 2 1 P3e{j}],-Ee{i}',1);
            lmiterm([LMI_obs2 2 2 -P4e{j}],1,-Ee{i},'s');
            
            lmiterm([LMI_obs2 1 1 -P3e{i}],1,Ae{j,k},'s');
            lmiterm([LMI_obs2 1 1 M1e{i}],1,-Ce{j,k},'s');
            lmiterm([LMI_obs2 1 1 P1e],1,2*beta);
            lmiterm([LMI_obs2 2 1 -P4e{i}],1,Ae{j,k});
            lmiterm([LMI_obs2 2 1 M2e{i}],1,-Ce{j,k});
            lmiterm([LMI_obs2 2 1 P1e],1,1);
            lmiterm([LMI_obs2 2 1 P3e{i}],-Ee{j}',1);
            lmiterm([LMI_obs2 2 2 -P4e{i}],1,-Ee{j},'s');
            
        end
    end
end

ejemplo=getlmis;
[tmin, xfeas]=feasp(ejemplo); %%manda a llamar el solver
factible=0;
if tmin<0
    factible=1;
    disp('FACTIBLE');
    P1=dec2mat(ejemplo,xfeas,P1);
    for j=1:r
       M{j}=dec2mat(ejemplo,xfeas,M{j}); 
       K{j}=M{j}*inv(P1);
    end
    P1e=dec2mat(ejemplo,xfeas,P1e);
    for j=1:re
        M1e{j}=dec2mat(ejemplo,xfeas,M1e{j});
        M2e{j}=dec2mat(ejemplo,xfeas,M2e{j});
        P3e{j}=dec2mat(ejemplo,xfeas,P3e{j});
        P4e{j}=dec2mat(ejemplo,xfeas,P4e{j});
    end
else
    factible=0;
end


