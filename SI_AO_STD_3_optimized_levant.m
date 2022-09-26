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

function [u,calculation_time] = SI_AO_STD_3_optimized_levant(f_k,h,L,reset_flag,order)
% Arbitrary differentiator order 3
%   Version 29-11-2019  Ok

persistent z0 z1 z2 z3 CT

    lambda0=3;
    lambda1=4.16;
    lambda2=3.06;
    lambda3=1.1;

    
    if reset_flag==1
        z0=0;
        z1=0;
        z2=0;
        z3=0;
        CT=0;
    end
    tic
    
%     z0_1=z0;
%     z1_1=z1;
%     z2_1=z2;
%     z3_1=z3;

    z3=h*(-lambda3*L*sign(z0-f_k))+z3;
    z2=h*(-lambda2*L^(3/4)*abs(z0-f_k)^(1/4)*sign(z0-f_k)+z3)+z2;
    z1=h*(-lambda1*L^(2/4)*abs(z0-f_k)^(2/4)*sign(z0-f_k)+z2)+z1+h^2*z3/2;
    z0=h*(-lambda0*L^(1/4)*abs(z0-f_k)^(3/4)*sign(z0-f_k)+z1)+z0+h^2*z2/2+h^3*z3/6;

    CT=CT+toc;
    
    if order==1
        u=z1;
    elseif order==2
        u=z2;
    elseif order==3
        u=z3;
    else
        disp('Arbitrary-order is not defined for n>3')
    end
    
	calculation_time=CT;    
end

