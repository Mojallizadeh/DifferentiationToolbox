% This file is part of the differentiation toolbox
% developed by Mohammad Rasool Mojallizadeh at INRIA, France
% Email: mohammad-rasool.mojallizadeh@inria.fr
%
% If you use this toolbox please cite the following report:
%
% Mohammad Rasool Mojallizadeh, Bernard Brogliato, Vincent Acary. Discrete-time differentiators:
% design and comparative analysis. 2020. hal-02960923
%
% This toolbox allows studying the behavior of several types of differentiator
% under different realistic conditions.
% 
% Copyright 2020 INRIA.
% 
% This program (the differentiation toolbox) is free software:
% you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%

clc;clear;close all;

o1=10*(rand(20,1)-0.5);
o2=10*(rand(20,1)-0.5);

subplot(2,2,1)

mu=0;

F_i=100;
alpha_i=5;

for j=1:20


h=1e-3;
lambda1=2;
lambda2=2;
reset_flag=1;


for k=1:100000
    
    f(k)=sin(h*k);
    df(k)=cos(h*k);
%     f(k)=0;
%     df(k)=0;
    
%     [X1(k),X2(k)] = STA_implicit_test(lambda1,lambda2,h,f(k),reset_flag,o1,o2);
    [X1(k),X2(k)] = high_degree_i_test(lambda1,lambda2,mu,h,f(k),reset_flag,o1(j),o2(j));
%     [X1(k),X2(k)] = quadratic_implicit_test(f(k),F_i,alpha_i,h,reset_flag,o1,o2);
    reset_flag=0;
    
    
end


plot(X1-f,X2-df)
hold on

plot(X1(1)-f(1),X2(1)-df(1),'x')

end

grid on


xlabel("$x_1-f$",'Interpreter','latex','FontSize',20)
ylabel("$x_2-df/dt$",'Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)
title("$\mu=0$",'Interpreter','latex','FontSize',20)


subplot(2,2,2)

mu=0.5;

F_i=100;
alpha_i=5;

for j=1:20


h=1e-3;
lambda1=2;
lambda2=2;
reset_flag=1;


for k=1:100000
    
    f(k)=sin(h*k);
    df(k)=cos(h*k);
%     f(k)=0;
%     df(k)=0;
    
%     [X1(k),X2(k)] = STA_implicit_test(lambda1,lambda2,h,f(k),reset_flag,o1,o2);
    [X1(k),X2(k)] = high_degree_i_test(lambda1,lambda2,mu,h,f(k),reset_flag,o1(j),o2(j));
%     [X1(k),X2(k)] = quadratic_implicit_test(f(k),F_i,alpha_i,h,reset_flag,o1,o2);
    reset_flag=0;
    
    
end


plot(X1-f,X2-df)
hold on

plot(X1(1)-f(1),X2(1)-df(1),'x')

end

grid on


xlabel("$x_1-f$",'Interpreter','latex','FontSize',20)
ylabel("$x_2-df/dt$",'Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)
title("$\mu=0.5$",'Interpreter','latex','FontSize',20)

subplot(2,2,3)

mu=1;

F_i=100;
alpha_i=5;

for j=1:20


h=1e-3;
lambda1=2;
lambda2=2;
reset_flag=1;


for k=1:100000
    
    f(k)=sin(h*k);
    df(k)=cos(h*k);
%     f(k)=0;
%     df(k)=0;
    
%     [X1(k),X2(k)] = STA_implicit_test(lambda1,lambda2,h,f(k),reset_flag,o1,o2);
    [X1(k),X2(k)] = high_degree_i_test(lambda1,lambda2,mu,h,f(k),reset_flag,o1(j),o2(j));
%     [X1(k),X2(k)] = quadratic_implicit_test(f(k),F_i,alpha_i,h,reset_flag,o1,o2);
    reset_flag=0;
    
    
end


plot(X1-f,X2-df)
hold on

plot(X1(1)-f(1),X2(1)-df(1),'x')

end

grid on


xlabel("$x_1-f$",'Interpreter','latex','FontSize',20)
ylabel("$x_2-df/dt$",'Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)
title("$\mu=1$",'Interpreter','latex','FontSize',20)

subplot(2,2,4)

mu=10;

F_i=100;
alpha_i=5;

for j=1:20


h=1e-3;
lambda1=2;
lambda2=2;
reset_flag=1;


for k=1:100000
    
    f(k)=sin(h*k);
    df(k)=cos(h*k);
%     f(k)=0;
%     df(k)=0;
    
%     [X1(k),X2(k)] = STA_implicit_test(lambda1,lambda2,h,f(k),reset_flag,o1,o2);
    [X1(k),X2(k)] = high_degree_i_test(lambda1,lambda2,mu,h,f(k),reset_flag,o1(j),o2(j));
%     [X1(k),X2(k)] = quadratic_implicit_test(f(k),F_i,alpha_i,h,reset_flag,o1,o2);
    reset_flag=0;
    
    
end


plot(X1-f,X2-df)
hold on

plot(X1(1)-f(1),X2(1)-df(1),'x')

end

grid on


xlabel("$x_1-f$",'Interpreter','latex','FontSize',20)
ylabel("$x_2-df/dt$",'Interpreter','latex','FontSize',20)
set(gca,'FontSize',20)
title("$\mu=10$",'Interpreter','latex','FontSize',20)