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

load data_Levant

% h=log10(h_vector);
h=h_vector;
% h=SNR;

% L_inf_STA_e=log10(L_inf_STA_e);
% L_inf_STA_i=log10(L_inf_STA_i);
% 
% L_inf_arbitrary_i_2=log10(L_inf_arbitrary_i_2);
% L_inf_arbitrary_i_3=log10(L_inf_arbitrary_i_3);
% L_inf_arbitrary_i_4=log10(L_inf_arbitrary_i_4);
% L_inf_arbitrary_i_5=log10(L_inf_arbitrary_i_5);
% L_inf_arbitrary_i_6=log10(L_inf_arbitrary_i_6);
% L_inf_arbitrary_i_7=log10(L_inf_arbitrary_i_7);
% 
% L_inf_arbitrary_e_2=log10(L_inf_arbitrary_e_2);
% L_inf_arbitrary_e_3=log10(L_inf_arbitrary_e_3);
% L_inf_arbitrary_e_4=log10(L_inf_arbitrary_e_4);
% L_inf_arbitrary_e_5=log10(L_inf_arbitrary_e_5);
% L_inf_arbitrary_e_6=log10(L_inf_arbitrary_e_6);
% L_inf_arbitrary_e_7=log10(L_inf_arbitrary_e_7);

% subplot(4,2,1)
% 
% 
% loglog(h,L_inf_arbitrary_e_1,'LineWidth',2)
% hold on
% loglog(h,L_inf_arbitrary_i_2,'LineWidth',2)
% legend("1st E-AO-STD","1st I-AO-STD",'FontSize',18)
% grid on
% set(gca,'FontSize',18)

subplot(3,2,1)

loglog(h,L_inf_arbitrary_e_2,'LineWidth',2)
hold on
loglog(h,L_inf_arbitrary_i_2,'LineWidth',2)
legend("2nd E-AO-STD","2nd I-AO-STD",'FontSize',18)
grid on
set(gca,'FontSize',18)
xlabel('$h$','Interpreter','latex','FontSize',18)

subplot(3,2,2)

loglog(h,L_inf_arbitrary_e_3,'LineWidth',2)
hold on
loglog(h,L_inf_arbitrary_i_3,'LineWidth',2)
legend("3rd E-AO-STD","3rd I-AO-STD",'FontSize',18)
grid on
set(gca,'FontSize',18)
xlabel('$h$','Interpreter','latex','FontSize',18)

subplot(3,2,3)

loglog(h,L_inf_arbitrary_e_4,'LineWidth',2)
hold on
loglog(h,L_inf_arbitrary_i_4,'LineWidth',2)
legend("4th E-AO-STD","4th I-AO-STD",'FontSize',18)
grid on
set(gca,'FontSize',18)

subplot(3,2,4)

loglog(h,L_inf_arbitrary_e_5,'LineWidth',2)
hold on
loglog(h,L_inf_arbitrary_i_5,'LineWidth',2)
legend("5th E-AO-STD","5th I-AO-STD",'FontSize',18)
grid on
set(gca,'FontSize',18)

subplot(3,2,5)

loglog(h,L_inf_arbitrary_e_6,'LineWidth',2)
hold on
loglog(h,L_inf_arbitrary_i_6,'LineWidth',2)
legend("6th E-AO-STD","6th I-AO-STD",'FontSize',18)
xlabel('$h$','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)

subplot(3,2,6)

loglog(h,L_inf_arbitrary_e_7,'LineWidth',2)
hold on
loglog(h,L_inf_arbitrary_i_7,'LineWidth',2)
legend("7th E-AO-STD","7th I-AO-STD",'FontSize',18)
xlabel('$h$','Interpreter','latex','FontSize',18)
grid on
set(gca,'FontSize',18)







