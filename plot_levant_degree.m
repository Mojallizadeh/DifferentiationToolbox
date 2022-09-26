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

clc;clear;close all;

Acc = [1e-12 1e-9 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1];
MNI = [23 23 22 21 20 18 25 32];
yyaxis left
semilogx(Acc,MNI,'LineWidth',2)

TNI = [3310 3183 2999 2750 1309 2493 3886 5227];
yyaxis right
semilogx(Acc,TNI,'LineWidth',2)
grid on

yyaxis left
% title('Plots with Different y-Scales','FontSize',18)
xlabel('Accuracy of the solver (Acc) (ms)','FontSize',18)
ylabel('MNI','FontSize',18)

yyaxis right
ylabel('TNI','FontSize',18)

set(gca,'FontSize',18)