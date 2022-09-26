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



function [u,calculation_time] = E_AO_STD_proper_6(f_k,h,L,reset_flag,order)
% Arbitrary differentiator order 3
%   Version 29-11-2019  Ok

persistent z0 z1 z2 z3 z4 z5 z6 CT

    lambda0=10;
    lambda1=47.69;
    lambda2=110.08;
    lambda3=101.96;
    lambda4=43.65;
    lambda5=9.91;
    lambda6=1.1;

    
    if reset_flag==1
        z0=0;
        z1=1;
        z2=0;
        z3=-1;
        z4=0;
        z5=1;
        z6=0;
        CT=0;
    end
    tic

    z0=h*(-lambda0*L^(1/7)*abs(z0-f_k)^(6/7)*sign(z0-f_k)+z1)+z0     +h^2*z2/2+h^3*z3/6+h^4*z4/24+h^5*z5/120+h^6*z6/720;   
    z1=h*(-lambda1*L^(2/7)*abs(z0-f_k)^(5/7)*sign(z0-f_k)+z2)+z1     +h^2*z3/2+h^3*z4/6+h^4*z5/24+h^5*z6/120; 
    z2=h*(-lambda2*L^(3/7)*abs(z0-f_k)^(4/7)*sign(z0-f_k)+z3)+z2     +h^2*z4/2+h^3*z5/6+h^4*z6/24; 
    z3=h*(-lambda3*L^(4/7)*abs(z0-f_k)^(3/7)*sign(z0-f_k)+z4)+z3     +h^2*z5/2+h^3*z6/6; 
    z4=h*(-lambda4*L^(5/7)*abs(z0-f_k)^(2/7)*sign(z0-f_k)+z5)+z4     +h^2*z6/2; 
    z5=h*(-lambda5*L^(6/7)*abs(z0-f_k)^(1/7)*sign(z0-f_k)+z6)+z5;
    z6=h*(-lambda6*L*sign(z0-f_k))+z6;
    

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

