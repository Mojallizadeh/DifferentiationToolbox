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

clear;clc

a=100;
bk=-100;

y=-500:0.001:500;

for j=1:length(y)
    
    if y(j)>bk
    f1(j)=(  (-a+sqrt(a^2-4*(bk-y(j))))/(2)    )^2;
    
    else
        
    f1(j)=-(  (a-sqrt(a^2+4*(bk-y(j))))/(2)    )^2;    
    end
    
end


plot(y,f1)
hold on
% plot(y,f2)
hold on