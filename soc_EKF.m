%soc(k)=soc(k-1)+T*i(k-1)/Q0;
%v(k)=3.272+2.837*soc(k)-8.452*soc(k)^2+7.012*soc(k)^3+13.06*soc(k)^4-24.11*soc(k)^5+10.62*soc(k)^6+T*R0*u(k)+Vc(k);
%C:h(k)=0.1+0.2*soc(k)^2
%% ----------------------完成----------------------
clc;clear;
soc(1)=0;
sochat(1)=0;
socne(1)=0;

Vc(1)=0;
Vchat(1)=0;
Vcne(1)=0;

y(1)=0.01;
p(1)=1;
psoc(1)=1;
psocne(1)=1;
pVc(1)=1;
pVcne(1)=1;

K(1)=0;
u=31;
Q=0.0001;
R=0.001;
T=0.01;
ESR=2.2*0.001;
R0=0.079;
R1=0.008038;
C1=33551.5256;
Q0=2.0962*3600;

A=[1, 0;
    0,1-T/C1/R1];
B=[T/Q0;T/C1];
D=R0;


for k=2:500
    soc(k)=soc(k-1)+T/Q0*u+normrnd(0,Q);
    Vc(k)=(1-T/C1/R1)*Vc(k-1)+T/C1*u+normrnd(0,Q);
    y(k) = 3.272+2.837*soc(k)-8.452*soc(k)^2+7.012*soc(k)^3+13.06*soc(k)^4-24.11*soc(k)^5+10.62*soc(k)^6+ T*R0*u+ Vc(k)+ normrnd(0,R);
    % H=[2.837-8.452*2*soc(k)+7.012*3*soc(k)^2+13.06*soc(k)^3*4-24.11*soc(k)^4*5+10.62*6*soc(k)^5, 1];
    % x(k)=A*x(k-1)+B*u+ones(2,1)*normrnd(0,Q);
    % y(k)=H*x(k) + D*u+normrnd(0,R);

    % xne(k)=xhat(k-1)+T*u;
    % yne(k)=xne(k)+T*ESR*u;
    socne(k)=sochat(k-1)+T/Q0*u;
    Vcne(k)=(1-T/C1/R1)*Vcne(k-1)+T/C1*u;
    yne(k)=f(socne(k))+Vcne(k)+T/C1*u;

    % pne(k)=1*p(k-1)*1+Q;
    % h(k)=1;
    psocne(k)=1*psoc(k-1)*1+Q;
    pVcne(k)=(1-T/C1/R1)*pVc(k-1)*(1-T/C1/R1)+Q;

    % K(k)=pne(k)*h(k)/(h(k)*pne(k)*h(k)+R);
    Ksoc(k)=psocne(k)*h_soc(socne(k))/(h_soc(socne(k))*psocne(k)*h_soc(socne(k))+R);
    KVc(k)=pVcne(k)*1/(1*pVcne(k)*1+R);

    % xhat(k)=xne(k)+K(k)*(y(k)-yne(k));
    sochat(k)=socne(k)+Ksoc(k)*(y(k)-yne(k));
    Vchat(k)=Vcne(k)+KVc(k)*(y(k)-yne(k));

    % p(k)=(1-K(k))*pne(k);
    psoc(k)=(1-Ksoc(k))*psocne(k);
    pVc(k)=(1-KVc(k))*pVcne(k);
end

t=1:500;
figure(1);
% plot(t(2:end),sochat(2:end),'r','LineWidth',2);
% legend("sochat");
% hold on;
% plot(t(2:end),Vchat(2:end),'g','LineWidth',2);
% legend("sochat","Vchat");
% hold on;
plot(t(2:end),soc(2:end),'r','LineWidth',2);
legend("sochat");
hold on;
plot(t(2:end),Vc(2:end),'g','LineWidth',2);
legend("sochat","Vchat");
hold on;

figure(2);
plot(t(2:end),y(2:end),'b','LineWidth',2);
title("输出");
legend("y");
hold on;    

function [result]=f(soc)
    result=3.272+2.837*soc-8.452*soc^2+7.012*soc^3+13.06*soc^4-24.11*soc^5+10.62*soc^6;
end

function [result]=h_soc(soc)
    result=2.837-8.452*2*soc+7.012*3*soc^2+13.06*soc^3*4-24.11*soc^4*5+10.62*6*soc^5;
end