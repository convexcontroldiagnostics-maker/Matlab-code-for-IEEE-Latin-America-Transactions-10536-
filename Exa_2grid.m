% Grid in Example 2, it compares the best results in Ricardo22

close all hidden;
clear all; clc;


% Parameters:
a=-0.5:0.05:1; b=-0.5:0.05:1;

Fig_sNQ=zeros(length(a),length(b));


for ii=1:length(a)
    for jj=1:length(b)
        [a(ii) b(jj)]
        
        
      
        %%% vertex matrices:
                    
        A{1}=[-1.15    0.1    1.8+b(jj);    0.3   -1.3   -0.5;   -0.1    0.8   -0.8];
        A{2}=[-1.2   -0.3   -0.1;    0.4   -0.6    0.3;   -0.2   -0.2   -0.2-a(ii)];
        B{1}=[0.6    1.2;    0.3    1.5-a(ii);   -0.6    1.3];
        B{2}=[-1.3    2.1;   -2.7    0.5;    1.5    1.6];
        C{1}=[ 0.4    1         0];
        C{2}=[ 0.8    1         0];
        E{1}=[1.05    0.7    0.7;   -0.1    1.1   -0.2;    0.1    0.5   0.9-a(ii)];
        E{2}=[0.9+b(jj)    0.8    0.77;   -0.9    1.1   -0.2;    0.4    0.5    0.6];
        
        r=2;

        [Feasp_sNQ]=Controller_Observer_Th1_Th2(A,B,C,E,r);
        if Feasp_sNQ==1
            Fig_sNQ(ii,jj)=1;
        end
        save dataGrid_Exa_Th1_Th2
    end
end

%

figure
%%%% Ploting results:
axes('FontSize',16,'FontName','Times New Roman');
xlabel('a','FontSize',16,'FontName','Times New Roman');
ylabel('b','FontSize',16,'FontName','Times New Roman');
box on;
hold on;
for i=1:length(a)
    for j=1:length(b)
        
   
        if Fig_sNQ(i,j)==1
            plot(a(i),b(j),'sb','LineWidth',1.6,'MarkerSize',12);
        end
   
    end
end

load('data_Ricardo.mat')

hold on;
for i=1:length(a)
    for j=1:length(b)
        
        if FigTh4_conf_3(i,j)==1
            plot(a(i),b(j),'xk','LineWidth',1.2,'MarkerSize',13);
        end
       if FigTh4_conf_9(i,j)==1
            plot(a(i),b(j),'xk','LineWidth',1.2,'MarkerSize',13);
        end
       
    end
end