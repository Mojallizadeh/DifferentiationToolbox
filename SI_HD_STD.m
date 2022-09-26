% This file is part of the differentiation toolbox
% developed by Mohammad Rasool Mojallizadeh at INRIA, France
% Email: mohammad-rasool.mojallizadeh@inria.fr
%
% If you use this toolbox please cite the following report:
%
% Mohammad Rasool Mojallizadeh, Bernard Brogliato, Vincent Acary. Discrete-time differentiators:
% design and comparative analysis. 2020. hal-02960923
%
% This toolbox allows studying the behavior of several types of
% differentiators under different realistic conditions.
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

function [u,calculation_time] = SI_HD_STD(L,mu,h,f_k,reset_flag,order,x0)
% Control signal of STA differentiato in explicit way
% Ver 24-11-2019  Ok

    lambda0=sqrt(L)*1.5;
    lambda1=L*1.1;

    persistent z0 z1 CT
    
    if reset_flag==1
        z0=x0(1);
        z1=x0(2);
        CT=0;
    end
    tic
    
     s0=(z0-f_k);
     
    alpha=(3/2)*h^2*lambda1*mu^2*abs(s0)^2+h*lambda0*mu*abs(s0)^(3/2)+h*lambda0*abs(s0)^(1/2)+(1+2*h^2*lambda1*mu)*abs(s0)+(1/2)*h^2*lambda1;
    
    psi0=(h*lambda0*mu*abs(s0)^(3/2)+2*abs(s0)+h*lambda0*abs(s0)^(1/2))/(alpha);
    
    psi1=(h*lambda0*mu*abs(s0)^(3/2)+abs(s0)+h*lambda0*abs(s0)^(1/2))/(alpha);
    
     z0=(-1+psi0)*s0+h*z1+f_k;
     z1=(-1+psi1)*s0/(h)+z1;
    
    
    
    if order==1
        u=z1;
    else
        u=0.1;
        disp("SI_HD_STD is not designed for higher-order differentiation")
    end
    CT=CT+toc;
    calculation_time=CT;
end
