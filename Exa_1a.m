function salida=Exa_1a(params)

global P1 P3 P4 K z_bound P1e P3e P4e M1e M2e tz_bound tzeta_bound
x1=params(1); x2=params(2); x3=params(3); x4=params(4);
x1g=params(5); x2g=params(6); x3g=params(7); x4g=params(8); tiempo=params(9);



%%% Controller scheduling functions
z(1)=x1^2;
z(2)=sinc(x3g);
w1(1)=(z_bound(1,2)-z(1))/(z_bound(1,2)-z_bound(1,1));  w1(2)=1-w1(1);
w2(1)=(z_bound(2,2)-z(2))/(z_bound(2,2)-z_bound(2,1));  w2(2)=1-w2(1);
W{1}=w1(1)*w2(1);
W{2}=w1(1)*w2(2);
W{3}=w1(2)*w2(1);
W{4}=w1(2)*w2(2);


Kw=0;
for j=1:numel(W)
    Kw=Kw+W{j}*K{j};
end

%%% Observer scheduling fucntions
tz(1)=x1^2;


tw1(1)=(tz_bound(1,2)-tz(1))/(tz_bound(1,2)-tz_bound(1,1));  tw1(2)=1-tw1(1);


tW{1}=tw1(1);
tW{2}=tw1(2);
M1ew=0; M2ew=0; P3ew=0; P4ew=0;
for j=1:numel(tW)
    P3ew=P3ew+tW{j}*P3e{j};
    P4ew=P4ew+tW{j}*P4e{j};
    M1ew=M1ew+tW{j}*M1e{j};
    M2ew=M2ew+tW{j}*M2e{j};
end

Pb=[P1e zeros(4); P3ew P4ew];
Mb=[M1ew; M2ew];
Lxgy=[eye(4) eye(4)]*inv(Pb)'*Mb;

u=Kw*[x1;x2g;x3g;x4g]+0*sin(tiempo);

%%% Plant
y=[x2+(x1^2+1)*x4;x1];
fx=[x1+x2+sin(x3)-0.1*x4; x1-2*x2; x1+x1^2*x2-0.3*x3;sin(x3)-x4];
gx=[x1^2+1;0;0;0];
dx=fx+gx*u;

%%% Observer:
yg=[x2g+(x1g^2+1)*x4g;x1g];
fxg=[x1+x2g+sin(x3g)-0.1*x4g; x1-2*x2g; x1+x1^2*x2g-0.3*x3g;sin(x3g)-x4g];
gxg=[x1^2+1;0;0;0];
dxg=fxg+gxg*u+Lxgy*(y-yg);

%error
e=[x1;x2;x3;x4]-[x1g;x2g;x3g;x4g];

V=[e(1) e(2) e(3) e(4)]*P1e*[e(1);e(2);e(3);e(4)];
salida=[dx;dxg;e;u];