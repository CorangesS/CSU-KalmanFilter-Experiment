%% ----------------------不知道能不能运行----------------------
R=2.2*0.001;
I=31*ones(1,100);
T=0.01;
t=1:100;
x_hat=zeros(1,100);
y_hat=zeros(1,100);
x    =zeros(1,100);
y    =zeros(1,100);

A=0;
B=1/310;
C=1;
D=R;
L=-1;

x_hat(1)=0;
y_hat(1)=0;
x(1)=0;
y(1)=0;

for i = 2:length(t)
    x_hat(i) = x_hat(i-1) + T*B*I(i);
    y_hat(i) = x_hat(i) + D*I(i);
    x(i) = x_hat(i-1) + T*B*I(i) + T*rand(1,1)*0.01;
    y(i) = x_hat(i-1) + T*B*I(i) + T*rand(1,1)*0.01;
    x_hat(i) = x_hat(i-1) + T*B*I(i) + L*(y(i)-y_hat(i));
    y_hat(i) = x_hat(i) + D*I(i);
end

figure;
plot(t,y);
