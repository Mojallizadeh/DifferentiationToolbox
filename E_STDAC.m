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


function [u,calculation_time,gamma1] = E_STDAC(alpha,eps,h,f_k,reset_flag,order,x0)
% Control signal of STA differentiato in explicit way
% Ver 24-11-2019  Ok

    lambda0=1.5;
    lambda1=1.1;
    persistent z0 z1 gamma CT
    
    if reset_flag==1
        z0=x0(1);
        z1=x0(2);
        gamma=1;
        CT=0;
    end
    
    tic

    s=(z0-f_k);
    
    if abs(s)>1
        gamma=h*(gamma/2)*alpha*abs(s)^(-1/2)+gamma;
    elseif abs(s)<1*eps
        gamma=h*((1/gamma)-1)+gamma;   
    else
        gamma=h*(gamma/2)*alpha*abs(s)+gamma;      
    end
    
	z0=-h*lambda0*gamma*abs(s)^(1/2)*sign(s)+h*z1+z0;
    z1=-h*lambda1*gamma^2*sign(s)+z1;
    
    if order==1
        u=z1;
    else
        u=0.1;
        disp("E_STDAC is not designed for higher-orders")
    end
    
    gamma1=gamma;
    CT=CT+toc;
    calculation_time=CT;
end
