%% -----------------------UKF模板-----------------------------------
clc;
clear;
%测量值模拟
T=0.01;%滤波周期
T_go=1;%滤波时间
N=T_go/T;%观测次数
t=0:T:T_go-T;%假定输出序列（供画图用）
x=zeros(2,N);
z=zeros(2,N);
x(:,1)=[0;0];%真值初始值
mu=[0;0];
Q=[0.0001,0;
    0,0.0001];
R=[0.0001,0;
    0,0.0001];
R0=0.079;
R1=0.008038;
C1=33551.5256;
Q0=2.0962*3600;


rng(1);
w=mvnrnd(mu,Q,N)';
v=mvnrnd(mu,R,N)';
for k=1:N-1
    x1=x(1,k);
    x2=x(2,k);
    w1=w(1,k);
    w2=w(2,k);
    % x(:,k+1)=[x1+T*x2+w1; -10*T*sin(x1)+(1-T)*x2+w2];
    x(:,k+1)=[(1-T/C1/R1)*x1+T/C1*31; x2+T/Q0*31];

    x1=x(1,k+1);
    x2=x(2,k+1);
    v1=v(1,k+1);
    v2=v(2,k+1);
    % z(:,k+1)=[2*sin(x1/2)+v1; x1/2+v2];
    z(:,k+1)=[x1+f(x2)+R0*31; 0];
end

%UKF估计
alpha=0.1;
beta=2;
kappa=1;
n=2;%状态维数
lamda=alpha^2*(n+kappa)-n;
xk=zeros(2,N);
zk=zeros(2,N);

%1.初始化
xk(:,1)=[0;0];
zk(:,1)=[0;0];
X0=mean(x);
Pk=mean((x-X0)*(x-X0)');

wm=zeros(1,5);
wc=zeros(1,5);
wm(1)=lamda/(n+lamda);
wc(1)=lamda/(n+lamda)+1-alpha^2+beta;
for i=2:2*n+1
    wm(i)=1/(2*(n+lamda));
    wc(i)=1/(2*(n+lamda));
end

for k=1:N-1
%2.计算k-1时刻采样点和权值
    xsigma=zeros(2,5);
    xsigma(:,1)=xk(:,k);
    xsigma(:,2)=xk(:,k)+sqrt(n+lamda)*sqrt(Pk(:,1));
    xsigma(:,3)=xk(:,k)+sqrt(n+lamda)*sqrt(Pk(:,2));
    xsigma(:,4)=xk(:,k)-sqrt(n+lamda)*sqrt(Pk(:,1));
    xsigma(:,5)=xk(:,k)-sqrt(n+lamda)*sqrt(Pk(:,2));
%3.计算k时刻的一步预测模型值
    xsigmapre=zeros(2,5);
    for i=1:5
        xsigmapre(:,i)=[(1-T/C1/R1)*xsigma(1,i)+T/C1*31; xsigma(2,i)+T/Q0*31];
    end
    xkpre=xsigmapre*wm';
    Pkpre=Q;
    for i=1:5
        Pkpre=Pkpre+wc(i)*(xsigmapre(:,i)-xkpre)*(xsigmapre(:,i)-xkpre)';
    end
%4.计算k时刻的一步预测样本点（重新采样）    
%     xsigma(:,1)=xkpre;
%     xsigma(:,2)=xkpre+sqrt(n+lamda)*sqrt(Pkpre(1));
%     xsigma(:,3)=xkpre+sqrt(n+lamda)*sqrt(Pkpre(2));
%     xsigma(:,4)=xkpre-sqrt(n+lamda)*sqrt(Pkpre(1));
%     xsigma(:,5)=xkpre-sqrt(n+lamda)*sqrt(Pkpre(2));
%5.计算k时刻的预测量测值
    zsigmapre=zeros(2,5);
    for i=1:5
       % zsigmapre(:,i)=[2*sin(xsigma(1,i)/2);0] ; 
       zsigmapre(:,i)=[xsigma(1,i)+f(xsigma(2,i))+R0*31;0];
    end
    zkpre=zsigmapre*wm';
    Pzzk=R;
    Pxzk=zeros(2);
    for i=1:5
        Pzzk=Pzzk+wc(i)*(zsigmapre(:,i)-zkpre)*(zsigmapre(:,i)-zkpre)';
        Pxzk=Pxzk+wc(i)*(xsigma(:,i)-xkpre)*(zsigmapre(:,i)-zkpre)';
    end
%6.量测更新
    Kk=Pxzk/Pzzk;
    xkpre=xkpre+Kk*(z(:,k+1)-zkpre);
    Pkpre=Pkpre-Kk*Pzzk*Kk';
    
    xk(:,k+1)=xkpre;
    zk(:,k+1)=zkpre;
    Pk=Pkpre;
end

figure;
% subplot(2,1,1);
plot(t,x(1,:),'-r',t,x(2,:),'-g');
% plot(t,x(2,:),'-g');
title("理想");
xlabel('t');
hold on;
legend('Vc理想','Soc理想');

figure;
% subplot(2,1,2);
plot(t,xk(1,:),'-r',t,xk(2,:),'-g');
% plot(t,xk(2,:),'-r');
title("UKF");
xlabel('t');
hold on;
legend('Vc','UKF');
% legend('Vc','Soc','V');

function [result]=f(soc)
    result=3.272+2.837*soc-8.452*soc^2+7.012*soc^3+13.06*soc^4-24.11*soc^5+10.62*soc^6;
end
function [result]=H(soc)
    result=2.837-8.452*soc*2+7.012*soc^2*3+13.06*soc^3*4-24.11*soc^4*5+10.62*soc^5*6;
end