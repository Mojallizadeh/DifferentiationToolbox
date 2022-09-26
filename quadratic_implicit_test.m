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

function [X1,X2] = quadratic_implicit_test(f_k,F_i,alpha_i,h,reset_flag,o1,o2)
%   Quadratic differentiator
%   Detailed explanation goes here

    persistent x1_1 x2_1
    
    if reset_flag==1
        x1_1 = o1;
        x2_1 = o2;
        CT= 0;
    end
    tic
    
	x2_star = sign(x1_1-f_k)*(F_i*h-sqrt(F_i^2*h^2+2*F_i*abs(x1_1-f_k)));

    x2_1 = x2_1-gsat(gsat(-alpha_i*h*F_i,x2_1,-h*F_i),x2_1-x2_star,gsat(h*F_i,x2_1,alpha_i*h*F_i));
    
    x1_1=h*x2_1+x1_1;

    u=x2_1;
    
 
    
    
	X1=x1_1;
    X2=x2_1;

end

