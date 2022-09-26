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


function [u,calculation_time] = E_QD(f_k,F_e,alpha_e,h,reset_flag,order,x0)
% Quadratic differentiator
%   Detailed explanation goes here

    persistent x1_1 x1_2 x1_3 x2_1 x2_2 x2_3 CT
    
    if reset_flag==1
        x1_1 = x0(1);
        x1_2 = x0(2);
        x1_3 = x0(3);        
        x2_1 = x0(2);
        x2_2 = x0(3); 
        x2_3 = x0(4);
        CT = 0;
    end
    tic
    
    
	s=2*F_e*(x1_1-f_k)+abs(x2_1)*x2_1;
    
    if s*x2_1>0
       x2_1=h*(-alpha_e*F_e*sign(s))+x2_1; 
    else
        x2_1=h*(-F_e*sign(s))+x2_1; 
    end

    
    x1_1=h*x2_1+x1_1;

    u=x2_1;
    
    if order==2 || order==3
    
        s=2*F_e*(x1_2-u)+abs(x2_2)*x2_2;
    
        if s*x2_2>0
            x2_2=h*(-alpha_e*F_e*sign(s))+x2_2; 
        else
            x2_2=h*(-F_e*sign(s))+x2_2; 
        end
    
    
        x1_2=h*x2_2+x1_2;

        u=x2_2;   
    
    end
    
    if order==3
    
        s=2*F_e*(x1_3-u)+abs(x2_3)*x2_3;
    
        if s*x2_3>0
            x2_3=h*(-alpha_e*F_e*sign(s))+x2_3; 
        else
            x2_3=h*(-F_e*sign(s))+x2_3; 
        end
    
    
        x1_3=h*x2_3+x1_3;

        u=x2_3;   
    
    end
    
    CT=CT+toc;
	calculation_time=CT; 
end

