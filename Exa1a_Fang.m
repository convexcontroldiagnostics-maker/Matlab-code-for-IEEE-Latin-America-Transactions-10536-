function salida=Exa1a_Fang(params)

x1=params(1); x2=params(2); x3=params(3); x4=params(4);

a=1.4; b=0.7;

%%% Controller scheduling functions
w1(1)=x1^2/a;  w1(2)=1-w1(1);
w2(1)=(b*sin(x3)-x3*sin(b))/(x3*b-x3*sin(b));  w2(2)=1-w2(1);
W{1}=w1(1)*w2(1);
W{2}=w1(1)*w2(2);
W{3}=w1(2)*w2(1);
W{4}=w1(2)*w2(2);

%%% Gains from the paper Fang2006_A New LMI-Based Approach to Relaxed Quadratic
%Stabilization of T–S Fuzzy Control Systems
%Section IV Examples

K{1}=[3.7255, 31.6739];
K{2}=[3.8230 32.3984];
K{3}=[7.5220, 44.0395];
K{4}=[7.5611, 44.4574];
Kw=0;
for j=1:numel(W)
    Kw=Kw+W{j}*K{j};
end

y=[x2+(x1^2+1)*x4;x1];
u=-Kw*y;

fx=[x1+x2+sin(x3)-0.1*x4; x1-2*x2; x1+x1^2*x2-0.3*x3;sin(x3)-x4];
gx=[x1^2+1;0;0;0];
dx=fx+gx*u;

salida=[dx;u];