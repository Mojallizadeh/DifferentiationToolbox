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

function [u,calculation_time] = SI_STD_sch2(L,h,f_k,reset_flag,order)
% Control signal of STA differentiato in explicit way
% Ver 24-11-2019  Ok

    lambda0=1.5;
    lambda1=1.1;

    persistent z0 z1 CT
    
    if reset_flag==1
        z0=0;z1=1;CT=0;
    end
    tic
   
        bk=f_k  -z0-h^1*z1...
            +h*lambda0*L^((1)/(2))*sign(z0-f_k)*abs(z0-f_k)^((1)/(2));
        
        if bk<-h^(2)*lambda1*L
            z1=h*(-lambda1*L)+z1;
            z0=h*(-lambda0*L^(1/2)*abs(z0-f_k)^(1/2)*sign(z0-f_k)+z1)+z0;            
             
        elseif bk>h^(2)*lambda1*L
            z1=h*(lambda1*L)+z1;
            z0=h*(-lambda0*L^(1/2)*abs(z0-f_k)^(1/2)*sign(z0-f_k)+z1)+z0;  
        else
            z1=bk/h+z1;
            z0=h*(-lambda0*L^(1/2)*abs(z0-f_k)^(1/2)*sign(z0-f_k)+z1)+z0;     
               
        end    
     
        CT=CT+toc;
        
	if order>1
        disp("Error in E-AO-STD")
	end

        calculation_time=CT;
        u=z1;
end
