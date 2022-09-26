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


x=-1:0.0001:1;

bk=3.5;

y00= x;
y01= x+bk;
y0_1= x-bk;


plot(y00,x,'color','b','LineWidth',2)
hold on
plot(y01,x,'color','b','LineWidth',2)
hold on
plot(y0_1,x,'color','b','LineWidth',2)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y20= x+abs(x).^(1/2).*sign(x);
y21= x+abs(x).^(1/2).*sign(x)+bk;
y2_1= x+abs(x).^(1/2).*sign(x)-bk;


plot(y10,x,'color','r','LineWidth',2)
hold on
plot(y11,x,'color','r','LineWidth',2)
hold on
plot(y1_1,x,'color','r','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y20= x+abs(x).^(1/3).*sign(x)+abs(x).^(2/3).*sign(x);
y21= x+abs(x).^(1/3).*sign(x)+abs(x).^(2/3).*sign(x)+bk;
y2_1= x+abs(x).^(1/3).*sign(x)+abs(x).^(2/3).*sign(x)-bk;


plot(y20,x,'color','g','LineWidth',2)
hold on
plot(y21,x,'color','g','LineWidth',2)
hold on
plot(y2_1,x,'color','g','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y30= x+abs(x).^(1/4).*sign(x)+abs(x).^(2/4).*sign(x)+abs(x).^(3/4).*sign(x);
y31= x+abs(x).^(1/4).*sign(x)+abs(x).^(2/4).*sign(x)+abs(x).^(3/4).*sign(x)+bk;
y3_1= x+abs(x).^(1/4).*sign(x)+abs(x).^(2/4).*sign(x)+abs(x).^(3/4).*sign(x)-bk;


plot(y30,x,'color','m','LineWidth',2)
hold on
plot(y31,x,'color','m','LineWidth',2)
hold on
plot(y3_1,x,'color','m','LineWidth',2)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

y40= x+abs(x).^(1/5).*sign(x)+abs(x).^(2/5).*sign(x)+abs(x).^(3/5).*sign(x)+abs(x).^(4/5).*sign(x);
y41= x+abs(x).^(1/5).*sign(x)+abs(x).^(2/5).*sign(x)+abs(x).^(3/5).*sign(x)+abs(x).^(4/5).*sign(x)+bk;
y4_1= x+abs(x).^(1/5).*sign(x)+abs(x).^(2/5).*sign(x)+abs(x).^(3/5).*sign(x)+abs(x).^(4/5).*sign(x)-bk;


plot(y40,x,'color','c','LineWidth',2)
hold on
plot(y41,x,'color','c','LineWidth',2)
hold on
plot(y4_1,x,'color','c','LineWidth',2)




xlim([-4 4])
% ylim([-.5 .5])
grid on