% Generate sample paths for the SDE by fBm
% dX_t = phi*f(X(t))dt +sigma(s)*dB_t^H where B_t^H is fractional brownian motion process

clc;
clear;
format long;
% randn('state',100);
H = 0.9; 
M = 2^13; 
T = 0.001*M;
dt=T/M;
t=0:dt:T;
N = 30;
X = zeros(M+1,N);
phi = normrnd(5,1,[N,1]);
d = @(x)      1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
f = @(x)     1+x;
X(1,:)=0;   %initial point
for i =1:M
 for j = 1:N
     [W,t]= fbm1d(H,M,T);
     X(i+1,j) = X(i,j) + phi(j)*f(X(i,j))*dt + d(X(i,j))*(W(i+1)-W(i));   
 end
end
save X305 X W phi
plot(t,X,'b-')


