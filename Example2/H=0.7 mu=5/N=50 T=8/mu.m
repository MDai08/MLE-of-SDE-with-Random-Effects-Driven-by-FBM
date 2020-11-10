clear;
clc;
%% 调用数据
V = load('X508.mat'); 
X = V.X;
W = V.W;
%% 设置参数
T = 8;   dt=0.001;   t=0:dt:T;
M = T/dt;   N = 50;
H = 0.7;
kH = 2*H*gamma(3/2-H)*gamma(H+1/2);
lamdaH = (2*H*gamma(3-2*H)*gamma(H+1/2))./(gamma(3/2-H));
KH = @(t,s)   inv(kH)*(s.^(1/2-H))*((t-s).^(1/2-H));
OH = @(t)     inv(lamdaH)*(t.^(2-2*H));
sigma = @(t)  1;
b = @(x,t)    1+x;
%% 计算PH
PH = zeros(M+1,N);
PH(1,:) = 0;
for j = 1:N
    PH(2,j) = 0;
end
for i =2:M
    for j = 1:N
        PH(i+1,j) = 0;
        for k=2:i
            PH(i+1,j) = PH(i+1,j)+(inv(kH)*(t(k).^(1/2-H))*((t(i+1)-t(k)).^(1/2-H))).*b(X(k,j)).*inv(sigma(t(k))).*(t(k+1)-t(k));
        end
    end
end


%% 计算QH
QH = zeros(M+1,N);
QH(1:2,:) = 0;
for j =1:N
    for i = 3:M+1
        QH(i,j) = (PH(i,j)-PH(i-1,j))/(inv(lamdaH)*(t(i).^(2-2*H)-t(i-1).^(2-2*H)));
    end
end


%% 计算Zt
Z = zeros(M+1,N);
Z(1,:) = 0;
for j = 1:N
    Z(2,j) = 0;
end
for i =2:M
    for j = 1:N
        Z(i+1,j) = 0;
        for k=2:i
            Z(i+1,j) = Z(i+1,j)+(inv(kH)*(t(k).^(1/2-H))*((t(i+1)-t(k)).^(1/2-H))).*inv(sigma(t(k))).*(X(k+1,j)-X(k,j));
        end
    end
end


%% 计算Uj
Y = zeros(M,N);
Omega = 1;
for j=1:N
    for i=1:M
        Y(i,j) = QH(i,j).*(Z(i+1,j)-Z(i,j));
    end
end

%% 计算Vj
P = zeros(M,N);
for j=1:N
    for i=1:M
        P(i,j) = (QH(i,j).^2).*(inv(lamdaH)*(t(i+1).^(2-2*H)-t(i).^(2-2*H)));
    end
end

save fbmdata508 QH Z kH lamdaH Y P
plot(t,Z,'b-');
%% 计算均值,方差
syms u w;
Ui= sum(Y);
Vi= sum(P);
Omega = 1;
u1 = 5;
u0 = double(solve(sum(Ui./(1+Omega^2.*Vi))-u*sum(Vi./(1+Omega^2.*Vi)),u))
w0 = double(solve(sum((u1-Ui./Vi).^2.*(Vi.^2./(1+w.*Vi).^2))-sum(Vi./(1+w.*Vi)),w))



