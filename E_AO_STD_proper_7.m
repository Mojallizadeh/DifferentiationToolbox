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


function [u,calculation_time] = E_AO_STD_proper_7(f_k,h,L,reset_flag,order)
% Arbitrary differentiator order 3
%   Version 29-11-2019  Ok

persistent z0 z1 z2 z3 z4 z5 z6 z7 CT

    lambda0=12;
    lambda1=84.14;
    lambda2=281.37;
    lambda3=455.40;
    lambda4=295.74;
    lambda5=88.78;
    lambda6=14.13;
    lambda7=1.1;

    
    if reset_flag==1
        z0=0;
        z1=1;
        z2=0;
        z3=-1;
        z4=0;
        z5=1;
        z6=0;
        z7=-1;        
        CT=0;
    end
    tic
   
    z0=h*(-lambda0*L^(1/8)*abs(z0-f_k)^(7/8)*sign(z0-f_k)+z1)+z0     +h^2*z2/2+h^3*z3/6+h^4*z4/24+h^5*z5/120+h^6*z6/720+h^7*z6/5040;     
    z1=h*(-lambda1*L^(2/8)*abs(z0-f_k)^(6/8)*sign(z0-f_k)+z2)+z1     +h^2*z3/2+h^3*z4/6+h^4*z5/24+h^5*z6/120+h^7*z6/720; 
    z2=h*(-lambda2*L^(3/8)*abs(z0-f_k)^(5/8)*sign(z0-f_k)+z3)+z2     +h^2*z4/2+h^3*z5/6+h^4*z6/24+h^5*z7/120; 
    z3=h*(-lambda3*L^(4/8)*abs(z0-f_k)^(4/8)*sign(z0-f_k)+z4)+z3     +h^2*z5/2+h^3*z6/6+h^4*z7/24; 
    z4=h*(-lambda4*L^(5/8)*abs(z0-f_k)^(3/8)*sign(z0-f_k)+z5)+z4     +h^2*z6/2+h^3*z7/6; 
    z5=h*(-lambda5*L^(6/8)*abs(z0-f_k)^(2/8)*sign(z0-f_k)+z6)+z5     +h^2*z7/2; 
    z6=h*(-lambda6*L^(7/8)*abs(z0-f_k)^(1/8)*sign(z0-f_k)+z7)+z6;
    z7=h*(-lambda7*L*sign(z0-f_k))+z7;
    

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

