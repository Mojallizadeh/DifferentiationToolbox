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

function [u,calculation_time] = STA_explicit(f_k,h,L,reset_flag,order)
% Arbitrary differentiator order 3
%   Version 29-11-2019  Ok

persistent z0 z1 CT

    lambda0=1.5;
    lambda1=1.1;

    if reset_flag==1
        z0=0;
        z1=1;
        CT=0;
    end
    tic
    
    z0_1=z0;
    z1_1=z1;

    z0=h*(-lambda0*L^(1/3)*abs(z0_1-f_k)^(2/3)*sign(z0_1-f_k)+z1_1)+z0_1;    
    z1=h*(-lambda1*L*sign(z0_1-f_k))+z1_1;
    
    CT=CT+toc;
    
    if order==1
        u=z1;
    else
        disp('E-STD is not defined for n>1')
    end
    
	calculation_time=CT;    
end

