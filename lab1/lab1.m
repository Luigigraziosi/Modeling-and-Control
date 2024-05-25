% Laboratorio 1 - Control of cyberphysical system
clear all
close all
clc

%Suppose that we dont' know the transfer function of the system,
% we had it only to be used to generate the plant input-output data by means of a simulation
% So with the functions rand and lsim we can generate input-output data,
% then we can compute the algorithm of System identification computing A
% and b (but not in our case) to raise the correct solution of the transfer
% function

%% Parameters
s = tf('s');
z = tf('z');

Gp = 100/(s^2+1.2*s+1);
Ts = 1;
Gd = c2d(Gp,Ts,'zoh');
n=2;

t = [7 10 25 30 45 1000 2000 5000 10000];
% t = [10];
%% System Identification Algorithm without noise
theta = zeros(5, length(t));
x = 1;
for H = t
    
    u = rand(H,1);
    y = lsim(Gd,u);
    
    A = zeros(H-2,2*n+1); 
    
    for i = 2: H-1
        A(i-1,:) = [ y(i) y(i-1) u(i+1) u(i) u(i-1)];
    end
    A;
    theta(:,x) = pinv(A)*y(3:H);
    x = x+1;
end

%% System Identification Algorithm with noise
theta1 = zeros(5,  length(t));
x = 1;
for H = t
    e = randn(H,1);
    [num, den] = tfdata(Gd, 'v');
    
    Gd_err = Gd + e/(den(1)*z^2+den(2)*z+den(3));
    Gd_err = Gd_err(1);
    u = rand(H,1);
    y = lsim(Gd_err,u);
    
    A = zeros(H-2,2*n+1); 
    
    for i = 2: H-1
        A(i-1,:) = [-y(i) -y(i-1) u(i+1) u(i) u(i-1)];
    end
    A;
    theta1(:,x) = pinv(A)*y(3:H);
    x=x+1;
end

%% System Identification Algorithm with output noise
theta2 = zeros(5, length(t));
x = 1;
for H = t
    u = rand(H,1);
    nu = randn(H,1);
    y = lsim(Gd,u) + nu;
    
    A = zeros(H-2,2*n+1); 
    
    for i = 2: H-1
        A(i-1,:) = [y(i) y(i-1) u(i+1) u(i) u(i-1)];
    end
    A;

    theta2(:,x) = pinv(A)*y(3:H);
    x = x+1;
end

theta
theta1
theta2








