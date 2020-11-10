% Generate sample paths for the SDE by fBm
% dX_t = f(X_t)dt +d1*dB_t^H where B_t^H is fractional brownian motion process

clc;
clear;
format long;
% randn('state',100)
H = 0.7; 
M = 2^13; 
T = 0.001*M;
dt=T/M;
t=0:dt:T;
N = 50;
X = zeros(M+1,N);
phi = normrnd(1,1,[N,1]);
d = @(x)      1;                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
X(1,:)=0;   %initial point
for i =1:M
 for j = 1:N
     [W,t]= fbm1d(H,M,T);
     X(i+1,j) = X(i,j) + phi(j)*dt + d(X(i,j))*(W(i+1)-W(i));   
 end
end
save X505 X W phi
plot(t,X,'b-')

