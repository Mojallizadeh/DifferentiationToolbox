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


function [u,calculation_time] = E_FDFF(omegas,omegaf,gama,f_k,h,reset_flag,order,x0)
%   Quadratic differentiator
%   Detailed explanation goes here

    persistent z1 z2 CT
    
    a1=gama;a2=gama*omegas;k1=2*omegaf;k2=omegaf^2;

    if reset_flag==1
        z1=x0(1);
        z2=x0(2);
        CT= 0;
    end
    tic
    
    z1=h*(z2-k1*(z1-f_k)-a1*sign(z1-f_k))+z1;
    z2=h*(-k2*(z1-f_k)-a2*sign(z1-f_k))+z2;
    
    out=z2;
    
    if order==1
        u=out;        
    else
        u=0;
    end
	    
	CT=CT+toc;
    
	calculation_time=CT;    

end

