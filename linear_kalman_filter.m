%% ----------------------不知道能不能运行----------------------
A=0;
B=1/310;
C=1;
Q=0.01;
R=0.1;
x_est=0;
P_est=0;
u=31;
ESR=2.2*0.001;
D=ESR;
y=C*x_est+D*u+rand(1,1)*0.01+rand(1,1)*0.01;

t=1:100;
result=zeros(1, 100);
result(1)=C*x_est+D*u;
for i = 2:length(t)
    [x_est,P_est] = t_linear_kalman_filter(A, B, C, Q, R, x_est, P_est, u, y);
    y=C*x_est+D*u+rand(1,1)*0.01+rand(1,1)*0.01;
    result(i)=C*x_est+D*u;
end

figure;
plot(t,result);

function [x_est, P_est] = t_linear_kalman_filter(A, B, C, Q, R, x_est, P_est, u, y)
    % A: 状态转移矩阵
    % B: 控制输入矩阵
    % C: 观测矩阵
    % Q: 过程噪声协方差
    % R: 观测噪声协方差
    % x_est: 初始状态估计
    % P_est: 初始误差协方差
    % u: 控制输入
    % y: 观测值

    % 预测步骤
    x_pred = A * x_est + B * u;
    P_pred = A * P_est * A' + Q;

    % 更新步骤
    K = P_pred * C' / (C * P_pred * C' + R); % 卡尔曼增益
    x_est = x_pred + K * (y - C * x_pred);
    P_est = (eye(size(K,1)) - K * C) * P_pred;

    % 返回更新后的状态估计和协方差估计
end
