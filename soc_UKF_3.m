clc;clear;
%UKF算法
t=0.01:0.01:1;
%生成真实的x与观测z
x=zeros(2,100);
z=zeros(2,100);
x(1,1)=0;
x(2,1)=0;
R0=0.079;
R1=0.008038;
C1=33551.5256;
Q0=2.0962*3600;
T=0.01;%滤波周期
%随便写一个非线性xk=f(xk-1) zk=h(xk)+R
for i=2:100
    % x(1,i)=sin(x(1,i-1))+2*cos(x(2,i-1));
    % x(2,i)=3*cos(x(1,i-1))+sin(x(2,i-1));
    x(1,i)=(1-T/C1/R1)*x(1,i-1)+T/C1*31;
    x(2,i)=x(2,i-1)+T/Q0*31;
    % z(1,i)=x(1,i)+x(2,i)^3+normrnd(0,1);
    % z(2,i)=x(1,i)^3+x(2,i)+normrnd(0,1);
    z(1,i)=x(1,i)+f(x(2,i))+R0*31;
    z(2,i)=0;
end

%设初值
X=zeros(2,100);
X(1,1)=0;
X(2,1)=0;
Pplus=[0.01,0;
    0,0.01];
Q=[0.0001,0;
    0,0.0001];
R=[0.01,0;
    0,0.01];


%设w(1,2n+1)
n=2;%n代表X的维数
w=zeros(n,2*n+1);
lamda = 2;
for i=1:2*n+1
    w(i)=1/(2*(n+lamda));
end
w(1)=lamda/(n+lamda);

%UKF算法实体
for i=2:100
    %拆x_sigma
    xsigma = zeros(n,2*n+1);
    L=chol(Pplus);
    %对正定矩阵进行分解
    xsigma(:,1)=X(:,i-1);
    for j=1:n
        xsigma(:,j+1)=xsigma(:,1)+sqrt(n+lamda)*L(:,j);
        xsigma(:,j+1+n)=xsigma(:,1)-sqrt(n+lamda)*L(:,j);
    end
    %预测步,生成xsigma_minus
    xsigmaminus=zeros(n,2*n+1);
    for j=1:2*n+1
        % xk=f(xk-1)
        % xsigmaminus(1,j)=sin(xsigma(1,j))+3*cos(xsigma(2,j));
        % xsigmaminus(2,j)=3*cos(xsigma(1,j))+sin(xsigma(2,j));
        xsigmaminus(1,j) = (1-T/C1/R1)*xsigma(1,j)+T/C1*31;
        xsigmaminus(2,j) = xsigma(2,j)+T/Q0*31;
    end
    %求期望协方差矩阵
    xhatminus=zeros(n,1);
    Pminus=zeros(n,n);
    for j=1:2*n+1
        xhatminus=xhatminus+w(j)*xsigmaminus(:,j);
    end
    for j=1:2*n+1
        Pminus=Pminus+w(j)*(xsigmaminus(:,j)-xhatminus)*(xsigmaminus(:,j)-xhatminus)';
    end
    Pminus=Pminus+Q;
    %更新步%再拆sigma点
    xsigma=zeros(n,2*n+1);
    xsigma(:,1)=xhatminus;
    L1=chol(Pminus);
    for j=1:n
        xsigma(:,j+1)=xsigma(:,1)+sqrt(n+lamda)*L1(:,j);
        xsigma(:,j+1+n)=xsigma(:,1)-sqrt(n+lamda)*L1(:,j);
    end
    %生成y,yhat
    yhat=zeros(n,1);
    for j=1:2*n+1
        % y(1,j)=xsigma(1,j)+xsigma(2,j)^3;
        % y(2,j)=xsigma(1,j)^3+xsigma(2,j);
        % yhat=yhat+w(j)*y(:,j);
        y(1,j)=xsigma(1,j)+f(xsigma(2,j))+R0*31;
        y(2,j)=0;
        yhat=yhat+w(j)*y(:,j);
    end
    %求Py Pxy
    Py=zeros(n,n);
    Pxy=zeros(n,n);
    for j=1:2*n+1
        Pxy=Pxy+w(j)*(xsigma(:,j)-xhatminus)*(y(:,j)-yhat)';
        Py=Py+w(j)*(y(:,j)-yhat)*(y(:,j)-yhat)';
    end
    Py=Py+R;
    %求卡尔曼增益
    K=Pxy*inv(Py);
    %观测数据
    Y=zeros(n,1);
    Y(1,1)=z(1,i);
    Y(2,1)=z(2,i);
    %更新
    X(:,i)=xhatminus+K*(Y-yhat);
    Pplus=Pminus+K*Py*K';
end

plot(t,x(2,:),'r',t,X(2,:),'b','LineWidth',3);
legend('真实','UKF估计');

function [result]=f(soc)
    result=3.272+2.837*soc-8.452*soc^2+7.012*soc^3+13.06*soc^4-24.11*soc^5+10.62*soc^6;
end
function [result]=H(soc)
    result=2.837-8.452*soc*2+7.012*soc^2*3+13.06*soc^3*4-24.11*soc^4*5+10.62*soc^5*6;
end