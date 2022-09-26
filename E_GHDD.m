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


function [u,calculation_time] = E_GHDD(f_k,h,L,reset_flag,order,x0)
% Arbitrary differentiator order 3
%   Version 5/17/2020  Ok

persistent z0 z1 z2 z3 CT

	lambda0=3;lambda1=4.16;lambda2=3.06;lambda3=1.1;   
    a12=-1/2;a13=1/3;a23=-1;
    
    if reset_flag==1
        z0=x0(1);
        z1=x0(2);
        z2=x0(3);
        z3=x0(4);
        CT=0;
    end
    tic

        
    z00=h*(-lambda0*L^(1/4)*abs(z0-f_k)^(3/4)*sign(z0-f_k)+z1)+z0 +h^2*z2/2 +h^3*z3/6;   
    z10=h*(-lambda1*L^(2/4)*abs(z0-f_k)^(2/4)*sign(z0-f_k)+z2)+z1 +h^2*z3/2              -a12*h^2*lambda2*L^(3/4)*abs(z0-f_k)^(1/4)*sign(z0-f_k)   -a13*h^3*L*lambda3*sign(z0-f_k) ; 
    z20=h*(-lambda2*L^(3/4)*abs(z0-f_k)^(1/4)*sign(z0-f_k)+z3)+z2                        -a23*h^2*lambda3*L*sign(z0-f_k)  ; 
    z30=h*(-lambda3*L*sign(z0-f_k))+z3; 
    
    z0=z00;
    z1=z10;
    z2=z20;
    z3=z30;
    

    CT=CT+toc;
    
    
    
    if order==0
        u=z0;
    elseif order==1
        u=z1;
    elseif order==2
        u=z2;
    elseif order==3
        u=z3;        
    else
        disp('Arbitrary-order is not defined for n>2')
    end
    
	calculation_time=CT;    
end

