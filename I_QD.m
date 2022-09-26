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


function [u,calculation_time] = I_QD(f_k,F_i,alpha_i,h,reset_flag,order,x0)
%   Quadratic differentiator
%   Detailed explanation goes here

    persistent x1_1 x1_2 x1_3 x2_1 x2_2 x2_3 CT
    
    if reset_flag==1
        x1_1 = x0(1);
        x1_2 = x0(2);
        x1_3 = x0(3);
        x2_1 = x0(2);
        x2_2 = x0(3);
        x2_3 = x0(4);
        CT= 0;
    end
    tic
    
	x2_star = sign(x1_1-f_k)*(F_i*h-sqrt(F_i^2*h^2+2*F_i*abs(x1_1-f_k)));

    x2_1 = x2_1-gsat(gsat(-alpha_i*h*F_i,x2_1,-h*F_i),x2_1-x2_star,gsat(h*F_i,x2_1,alpha_i*h*F_i));
    
    x1_1=h*x2_1+x1_1;

    u=x2_1;
    
 
    
    if order==2 || order==3
        
        x2_star = sign(x1_2-u)*(F_i*h-sqrt(F_i^2*h^2+2*F_i*abs(x1_2-u)));

        x2_2 = x2_2-gsat(gsat(-alpha_i*h*F_i,x2_2,-h*F_i),x2_2-x2_star,gsat(h*F_i,x2_2,alpha_i*h*F_i));

        x1_2=h*x2_2+x1_2;

        u=x2_2;    
        
    end
    
    if order==3
        
        x2_star = sign(x1_3-u)*(F_i*h-sqrt(F_i^2*h^2+2*F_i*abs(x1_3-u)));

        x2_3 = x2_3-gsat(gsat(-alpha_i*h*F_i,x2_3,-h*F_i),x2_3-x2_star,gsat(h*F_i,x2_3,alpha_i*h*F_i));

        x1_3=h*x2_3+x1_3;

        u=x2_3;             
        
    end

    
    CT=CT+toc;
    
	calculation_time=CT;    

end

