Vocv = [3.2976,3.4795,3.5715,3.6069,3.6344,3.6756,3.7657,3.8542,3.9544,4.0643,4.1784];
Soc = [0,10,20,30,40,50,60,70,80,90,100];

[slope, intercept] = Soc_Q_nihe_fun(Soc, Vocv);
fprintf('线性拟合函数为 y = %f * x + %f\n', slope, intercept);

function [slope, intercept] = Soc_Q_nihe_fun(x, y)
    % 检查输入向量的大小
    if length(x) ~= length(y) || size(x, 1) ~= 1 || size(y, 1) ~= 1
        error('x 和 y 必须是相同长度的一维向量');
    end

    % 使用 polyfit 进行线性拟合
    % polyfit 返回一次多项式的系数：y = slope * x + intercept
    p = polyfit(x, y, 1);

    % 提取斜率和截距
    slope = p(1);
    intercept = p(2);
end
