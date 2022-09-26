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


function [u,calculation_time] = I_FDFF(omegas,omegaf,gama,ro,f_k,h,reset_flag,order,x0)
%   Quadratic differentiator
%   Detailed explanation goes here

    persistent z1 z2 CT D B
    
    a1=gama;a2=gama*omegas;k1=2*omegaf;k2=omegaf^2;

    if reset_flag==1
        z1=x0(1);
        z2=x0(2);
        CT= 0;
        D=1+h*k1+h^2*k2;
        B=h*(a1+h*a2)/D;
    end
    tic
    
    z1m=z1;
    z2m=z2;
    
	e_star=(z1m+h*z2m-f_k)/D;
    s=psi_kikk(e_star/B,ro);
    
    e1=e_star-B*s;
    z1=f_k+e1;
    z2=(z1-z1m)/h+k1*e1+a1*s;
    
    
    out=z2;
    
    if order==1
        u=out;        
    else
        u=0;
    end
	    
	CT=CT+toc;
    
	calculation_time=CT;    

end

