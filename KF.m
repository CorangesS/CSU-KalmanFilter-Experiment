%% ----------------------模板----------------------
clc;
clear;

Z = 24 + randn(1, 4000);     %假设测量值是24加上高斯噪声N(0,1)
X = zeros(1, 4000);          %最优估的初始值，后面会迭代更新
X(1)=24;

P = zeros(1, 4000);          %预测模型的协方差，注意这里并不是一个对称矩阵
Q = 0.001;                   %相信观测转移矩阵，误差设置的较小
R = 0.01;                    %观测模型的协方差

for t=2:4000
    X_(t) = X(t-1);    
    P_(t) = P(t-1) + Q;
    K(t) = P_(t) / (P(t) + R);
    X(t) = X_(t) + K(t) * (Z(t) - X_(t));
    P(t) = (1 - K(t)) * P_(t);
end

 
figure;
plot(Z, 'k+');
hold on;

plot(X, 'b-', 'linewidth', 1);
hold on;
plot(24*ones(1, 4000), 'r-', 'linewidth', 1);