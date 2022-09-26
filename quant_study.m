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

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%

clc;clear;close all

h=1e-6;
q=10^-6;
t=0:h:10;
x=sin(t);
y=quant(x,q);
subplot(2,2,1)
plot(t,y,'LineWidth',2)
grid on
xlabel('Time ($h=10^{-6}$s)','Interpreter','latex','FontSize',18)
ylabel('Resolution = $10^{-6}$','Interpreter','latex','FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(A)','FontSize',18) 
set(gca,'FontSize',18)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=1e-6;
q=10^-1;
t=0:h:10;
x=sin(t);
y=quant(x,q);
subplot(2,2,2)
plot(t,y,'LineWidth',2)
grid on
xlabel('Time ($h=10^{-6}$s)','Interpreter','latex','FontSize',18)
ylabel('Resolution = $10^{-1}$','Interpreter','latex','FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(B)','FontSize',18) 
set(gca,'FontSize',18)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=0.1;
q=10^-6;
t=0:h:10;
x=sin(t);
y=quant(x,q);
subplot(2,2,3)
plot(t,y,'LineWidth',2)
grid on
xlabel('Time ($h=10^{-1}$s)','Interpreter','latex','FontSize',18)
ylabel('Resolution = $10^{-6}$','Interpreter','latex','FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(C)','FontSize',18) 
set(gca,'FontSize',18)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

h=0.1;
q=10^-1;
t=0:h:10;
x=sin(t);
y=quant(x,q);
subplot(2,2,4)
plot(t,y,'LineWidth',2)
grid on
xlabel('Time ($h=10^{-1}$s)','Interpreter','latex','FontSize',18)
ylabel('Resolution = $10^{-1}$','Interpreter','latex','FontSize',18)
ylim=get(gca,'ylim');
xlim=get(gca,'xlim');
text(xlim(end),ylim(end),'(D)','FontSize',18) 
set(gca,'FontSize',18)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%