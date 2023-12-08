    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  容积Kalman滤波
    %  状态方程：x(:,k+1) = F * x(:,k) + [T/C1,T/Q0]*1+ [sqrt(Q) * randn;sqrt(Q) * randn]; 
    %  观测方程：z(k+1) = f(x(1,k+1)) + x(2,k+1) + sqrt(R) * randn;
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 %% ----------------------完成----------------------
    clc; clear;
    n=501;
    tf = 500;                                     % 模拟长度 
    x=zeros(2,n);
    z=zeros(1,n);
    x(:,1) =[0;0];                              % 初始状态 

    x_ckf=zeros(2,n);
    % x_estimate(:,1) = [1;0.1];                  %状态的估计
    x_ckf(:,1)=[0;0];
    % e_x_estimate = x_estimate(:,1);             %EKF的初始估计
    xhat=x_ckf(:,1);
    yhat=zeros(1,n);
    x_e_error=zeros(1,n);
    x_c_error=zeros(1,n);
    z_e_error=zeros(1,n);
    z_c_error=zeros(1,n);

    u=31;
    Q=0.0001;
    R=0.01;
    T=0.01;
    ESR=2.2*0.001;
    R0=0.079;
    R1=0.008038;
    C1=33551.5256;
    Q0=2.0962*3600;

    P =[0.00001,0;
    0,0.00001];                        %初始估计方差
    Pplus=P;
    F=[1-T/C1/R1,0;
        0,1];
    Gamma=[1;1];
    w=0.25;  
    kesi=sqrt(2)*[1,0,-1,0;
        0,1,0,-1];

    for k = 1 : tf 
       % 模拟系统 
       x(:,k+1) = F * x(:,k) + [T/C1;T/Q0]*31 + Gamma * sqrt(Q) * randn;      %状态值 
       %x(:,k+1) = F * x(:,k) +[sqrt(Q) * randn;0]; 
       z(k+1) = f(x(1,k+1)) + x(2,k+1) + R0*31 + sqrt(R) * randn;      %观测值
    end
    
    for k = 1 : tf 
        %Cubature卡尔曼滤波器
        %%%%%（1）求协方差矩阵平方根
        S=chol(Pplus,'lower');
        %%%%%（2）计算求容积点
        rjpoint(:,1)=S*kesi(:,1)+xhat;
        rjpoint(:,2)=S*kesi(:,2)+xhat;
        rjpoint(:,3)=S*kesi(:,3)+xhat;
        rjpoint(:,4)=S*kesi(:,4)+xhat;
        %%%%%（3）传播求容积点
        Xminus(:,1)=F*rjpoint(:,1);                           %容积点经过非线性函数后的值
        Xminus(:,2)=F*rjpoint(:,2);
        Xminus(:,3)=F*rjpoint(:,3);
        Xminus(:,4)=F*rjpoint(:,4);
        %%%%（4）状态预测
        xminus=w*Xminus(:,1)+w*Xminus(:,2)+w*Xminus(:,3)+w*Xminus(:,4);
        %%%%(5)状态预测协方差阵
        Pminus=w*(Xminus(:,1)*Xminus(:,1)'+Xminus(:,2)*Xminus(:,2)'+Xminus(:,3)*Xminus(:,3)'+Xminus(:,4)*Xminus(:,4)')-xminus*xminus'+Gamma * Q* Gamma';
       %Pminus=w*(Xminus(:,1)*Xminus(:,1)'+Xminus(:,2)*Xminus(:,2)'+Xminus(:,3)*Xminus(:,3)'+Xminus(:,4)*Xminus(:,4)')-xminus*xminus'+[Q,0;0,0]; 
       %%%%观测更新
        %%%%%（1）矩阵分解
        Sminus=chol(Pminus,'lower');
        %%%%%（2）计算求容积点
        rjpoint1(:,1)=Sminus*kesi(:,1)+xminus;
        rjpoint1(:,2)=Sminus*kesi(:,2)+xminus;
        rjpoint1(:,3)=Sminus*kesi(:,3)+xminus;
        rjpoint1(:,4)=Sminus*kesi(:,4)+xminus;
        %%%%%（3）传播求容积点
        % Z(1)=atan(0.1*rjpoint1(1,1));
        % Z(2)=atan(0.1*rjpoint1(1,2));
        % Z(3)=atan(0.1*rjpoint1(1,3));
        % Z(4)=atan(0.1*rjpoint1(1,4));
        Z(1)=f(rjpoint1(1,1)) + rjpoint1(2,1) + R0*31;
        Z(2)=f(rjpoint1(1,2)) + rjpoint1(2,2) + R0*31;
        Z(3)=f(rjpoint1(1,3)) + rjpoint1(2,3) + R0*31;
        Z(4)=f(rjpoint1(1,4)) + rjpoint1(2,4) + R0*31;
        % z(k+1) = f(x(1,k+1)) + x(2,k+1) + R0*1 + sqrt(R) * randn;      %观测值
        %%%%%%%（4）观测预测
        zhat=w*(Z(1)+Z(2)+Z(3)+Z(4));
        %%%%(5)观测预测协方差阵
        %Pzminus=w*(Z(:,1)*Z(:,1)'+Z(:,2)*Z(:,2)'+Z(:,3)*Z(:,3)'+Z(:,4)*Z(:,4)')-zhat*zhat'+[R,0;0,Q];
        Pzminus=w*(Z(1)^2+Z(2)^2+Z(3)^2+Z(4)^2)-zhat^2+R;
        %%%%(6)互协方差阵
        Pxzminus=w*(rjpoint1(:,1)*Z(1)+rjpoint1(:,2)*Z(2)+rjpoint1(:,3)*Z(3)+rjpoint1(:,4)*Z(4))-xminus*zhat;
        %%%%(7)计算卡尔曼增益
        K=Pxzminus/Pzminus;
        %%%%(8)状态更新
        xhat=xminus+K*(z(k+1)-zhat);
        %%%%(9)状态协方差矩阵更新
        Pplus=Pminus-K*Pzminus*K';
        
        x_ckf(:,k+1)=xhat;
        yhat(k+1)=f(x_ckf(1,k+1)) + x_ckf(2,k+1) + R0*31;
    end
    
    
    t = 0 : tf;
    figure;
    plot(t,x(1,:),'k.',t,x_ckf(1,:),'g');
    legend('真实值','CKF估计值');
    figure;
    plot(t,yhat(1,:),'k.');

function [result]=f(soc)
    result=3.272+2.837*soc-8.452*soc^2+7.012*soc^3+13.06*soc^4-24.11*soc^5+10.62*soc^6;
end