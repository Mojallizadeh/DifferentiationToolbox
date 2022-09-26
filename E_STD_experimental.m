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


function [u,x_11] = E_STD_experimental(L,h,f_k,reset_flag)
% Control signal of STA differentiato in explicit way
% Ver 24-11-2019  Ok

    lambda1=sqrt(L)*1.5;
    lambda2=L*1.1;

    persistent v_1 x_1
    
    if reset_flag==1
        x_1 = 0;
        v_1 = 1;
    end
    

    
        u_1=v_1-lambda1*sqrt(abs(x_1-f_k))*sign(x_1-f_k);
        v_1=h*(-lambda2*sign(x_1-f_k))+v_1;
        x_1=h*u_1+x_1;
        
        u=u_1;
        x_11=x_1;

end
