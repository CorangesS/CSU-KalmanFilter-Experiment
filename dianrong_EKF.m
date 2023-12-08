%soc(k)=soc(k-1)+T*i(k-1)/Q0;
clear;
clc;
x(1)=0;
xhat(1)=0;
xne(1)=0;
y(1)=0.01;
p(1)=1;
K(1)=0;
u=31;
Q=0.001;
R=0.01;
T=0.01;
dianrong_c=310;
ESR=2.2*0.001;

for k=2:500
    x(k)=x(k-1)+T/dianrong_c*u+normrnd(0,Q);
    y(k)=x(k)+ T*ESR*u+normrnd(0,R);

    xne(k)=xhat(k-1)+T*u;
    yne(k)=xne(k)+T*ESR*u;
    
    pne(k)=1*p(k-1)*1+Q;
    h(k)=1;
    K(k)=pne(k)*h(k)/(h(k)*pne(k)*h(k)+R);
    xhat(k)=xne(k)+K(k)*(y(k)-yne(k));
    p(k)=(1-K(k))*pne(k);
end

t=1:500;
figure(1);
plot(t,xhat,'r','LineWidth',2);
legend("x_'预测");
hold on;
plot(t,x,'g','LineWidth',2);
legend("x_h_a_t","x状态");
hold on;

figure(2);
plot(t,y,'b','LineWidth',2);
title("输出");
legend("y");
hold on;    