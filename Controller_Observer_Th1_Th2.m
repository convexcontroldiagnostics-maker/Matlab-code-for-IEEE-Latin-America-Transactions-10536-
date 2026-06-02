%%% LMI conditions in Theorems 1 and 2 

function [factible]=Controller_Observer_Th1_Th2(A,B,C,E,r)


factible=0;


[n m]=size(B{1});
[q n]=size(C{1});

setlmis([]); %%%inicializa un grupo de LMIs

%decision variables
P1=lmivar(1,[n 1]);
P1e=lmivar(1,[n 1]);
for j=1:r
    P3{j}=lmivar(2,[n n]);
    P4{j}=lmivar(2,[n n]);
    M{j}=lmivar(2,[m n]);
    P3e{j}=lmivar(2,[n n]);
    P4e{j}=lmivar(2,[n n]);
    M1e{j}=lmivar(2,[n q]);
    M2e{j}=lmivar(2,[n q]);
end

%%% LMIs from Theorem 1
LMI0=newlmi;
lmiterm([-LMI0 1 1 P1],1,1);

for i=1:r
    for j=1:r
        LMI2=newlmi;
        lmiterm([LMI2 1 1 P3{i}],1,2/(r-1),'s');
        lmiterm([LMI2 2 1 P1],A{i},2/(r-1));
        lmiterm([LMI2 2 1 M{i}],B{i},2/(r-1));
        lmiterm([LMI2 2 1 P3{i}],-E{i},2/(r-1));
        lmiterm([LMI2 2 1 -P4{i}],1,2/(r-1));
        lmiterm([LMI2 2 2 P4{i}],-E{i},2/(r-1),'s');
        
        lmiterm([LMI2 1 1 P3{j}],1,1,'s');
        lmiterm([LMI2 2 1 P1],A{i},1);
        lmiterm([LMI2 2 1 M{j}],B{i},1);
        lmiterm([LMI2 2 1 P3{j}],-E{i},1);
        lmiterm([LMI2 2 1 -P4{j}],1,1);
        lmiterm([LMI2 2 2 P4{j}],-E{i},1,'s');
        
        lmiterm([LMI2 1 1 P3{i}],1,1,'s');
        lmiterm([LMI2 2 1 P1],A{j},1);
        lmiterm([LMI2 2 1 M{i}],B{j},1);
        lmiterm([LMI2 2 1 P3{i}],-E{j},1);
        lmiterm([LMI2 2 1 -P4{i}],1,1);
        lmiterm([LMI2 2 2 P4{i}],-E{j},1,'s');
        
    end
end


%%% LMIs from Theorem 2
LMI1=newlmi;
lmiterm([-LMI1 1 1 P1e],1,1);

for i=1:r
    for j=1:r
        LMI3=newlmi;
        lmiterm([LMI3 1 1 -P3e{i}],2/(r-1),A{i},'s');
        lmiterm([LMI3 1 1 M1e{i}],2/(r-1),-C{i},'s');
        lmiterm([LMI3 2 1 -P4e{i}],2/(r-1),A{i});
        lmiterm([LMI3 2 1 M2e{i}],2/(r-1),-C{i});
        lmiterm([LMI3 2 1 P1e],2/(r-1),1);
        lmiterm([LMI3 2 1 P3e{i}],-E{i}',2/(r-1));
        lmiterm([LMI3 2 2 -P4e{i}],2/(r-1),-E{i},'s');
                
        lmiterm([LMI3 1 1 -P3e{j}],1,A{i},'s');
        lmiterm([LMI3 1 1 M1e{j}],1,-C{i},'s');
        lmiterm([LMI3 2 1 -P4e{j}],1,A{i});
        lmiterm([LMI3 2 1 M2e{j}],1,-C{i});
        lmiterm([LMI3 2 1 P1e],1,1);
        lmiterm([LMI3 2 1 P3e{j}],-E{i}',1);
        lmiterm([LMI3 2 2 -P4e{j}],1,-E{i},'s');
        
        lmiterm([LMI3 1 1 -P3e{i}],1,A{j},'s');
        lmiterm([LMI3 1 1 M1e{i}],1,-C{j},'s');
        lmiterm([LMI3 2 1 -P4e{i}],1,A{j});
        lmiterm([LMI3 2 1 M2e{i}],1,-C{j});
        lmiterm([LMI3 2 1 P1e],1,1);
        lmiterm([LMI3 2 1 P3e{i}],-E{j}',1);
        lmiterm([LMI3 2 2 -P4e{i}],1,-E{j},'s');
        
    end
end


ejemplo=getlmis;
[tmin, xfeas]=feasp(ejemplo); %%manda a llamar el solver
factible=0;
if tmin<0
    factible=1;
    disp('FACTIBLE new');
    P1e=dec2mat(ejemplo,xfeas,P1e);
    for j=1:r
        P3e{j}=dec2mat(ejemplo,xfeas,P3e{j});
        P4e{j}=dec2mat(ejemplo,xfeas,P4e{j});
        M1e{j}=dec2mat(ejemplo,xfeas,M1e{j});
        M2e{j}=dec2mat(ejemplo,xfeas,M2e{j});
    end
else
    factible=0;
end


