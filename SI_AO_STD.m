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

function [u,calculation_time] = SI_AO_STD(f_k,h,L,reset_flag,output_order,differentiator_order)
% Arbitrary differentiator order 3
%   Version 29-11-2019  Ok

persistent z0 z1 z2 z3 z4 z5 z6 z7 CT

    if reset_flag==1
        z0=0;z1=1;z2=0;z3=-1;z4=0;z5=1;z6=0;z7=-1;CT=0;
    end
    tic
    
    z0_1=z0;z1_1=z1;z2_1=z2;z3_1=z3;z4_1=z4;z5_1=z5;z6_1=z6;z7_1=z7;
    
    if differentiator_order==7
        lambda0=12;lambda1=84.14;lambda2=281.37;lambda3=455.40;lambda4=295.74;lambda5=88.78;lambda6=14.13;lambda7=1.1;
        z0=h*(-lambda0*L^(1/8)*abs(z0_1-f_k)^(7/8)*sign(z0_1-f_k)+z1_1)+z0_1;     
        z1=h*(-lambda1*L^(2/8)*abs(z0_1-f_k)^(6/8)*sign(z0_1-f_k)+z2_1)+z1_1; 
        z2=h*(-lambda2*L^(3/8)*abs(z0_1-f_k)^(5/8)*sign(z0_1-f_k)+z3_1)+z2_1; 
        z3=h*(-lambda3*L^(4/8)*abs(z0_1-f_k)^(4/8)*sign(z0_1-f_k)+z4_1)+z3_1; 
        z4=h*(-lambda4*L^(5/8)*abs(z0_1-f_k)^(3/8)*sign(z0_1-f_k)+z5_1)+z4_1; 
        z5=h*(-lambda5*L^(6/8)*abs(z0_1-f_k)^(2/8)*sign(z0_1-f_k)+z6_1)+z5_1; 
        z6=h*(-lambda6*L^(7/8)*abs(z0_1-f_k)^(1/8)*sign(z0_1-f_k)+z7_1)+z6_1;
        z7=h*(-lambda7*L*sign(z0_1-f_k))+z7_1;
    
    elseif differentiator_order==6
        lambda0=10;lambda1=47.69;lambda2=110.08;lambda3=101.96;lambda4=43.65;lambda5=9.91;lambda6=1.1;
        z0=h*(-lambda0*L^(1/7)*abs(z0_1-f_k)^(6/7)*sign(z0_1-f_k)+z1_1)+z0_1;   
        z1=h*(-lambda1*L^(2/7)*abs(z0_1-f_k)^(5/7)*sign(z0_1-f_k)+z2_1)+z1_1; 
        z2=h*(-lambda2*L^(3/7)*abs(z0_1-f_k)^(4/7)*sign(z0_1-f_k)+z3_1)+z2_1; 
        z3=h*(-lambda3*L^(4/7)*abs(z0_1-f_k)^(3/7)*sign(z0_1-f_k)+z4_1)+z3_1; 
        z4=h*(-lambda4*L^(5/7)*abs(z0_1-f_k)^(2/7)*sign(z0_1-f_k)+z5_1)+z4_1; 
        z5=h*(-lambda5*L^(6/7)*abs(z0_1-f_k)^(1/7)*sign(z0_1-f_k)+z6_1)+z5_1;
        z6=h*(-lambda6*L*sign(z0_1-f_k))+z6_1;
        
    elseif differentiator_order==5
        lambda0=7;lambda1=23.72;lambda2=32.24;lambda3=20.26;lambda4=6.75;lambda5=1.1;        
        z0=h*(-lambda0*L^(1/6)*abs(z0_1-f_k)^(5/6)*sign(z0_1-f_k)+z1_1)+z0_1;    
        z1=h*(-lambda1*L^(2/6)*abs(z0_1-f_k)^(4/6)*sign(z0_1-f_k)+z2_1)+z1_1;
        z2=h*(-lambda2*L^(3/6)*abs(z0_1-f_k)^(3/6)*sign(z0_1-f_k)+z3_1)+z2_1;
        z3=h*(-lambda3*L^(4/6)*abs(z0_1-f_k)^(2/6)*sign(z0_1-f_k)+z4_1)+z3_1;
        z4=h*(-lambda4*L^(5/6)*abs(z0_1-f_k)^(1/6)*sign(z0_1-f_k)+z5_1)+z4_1;
        z5=h*(-lambda5*L*sign(z0_1-f_k))+z5_1;
    
	elseif differentiator_order==4
        lambda1=5;lambda2=10.03;lambda3=9.3;lambda4=4.57;lambda5=1.1;
        z0=h*(-lambda1*L^(1/5)*abs(z0_1-f_k)^(4/5)*sign(z0_1-f_k)+z1_1)+z0_1;
        z1=h*(-lambda2*L^(2/5)*abs(z0_1-f_k)^(3/5)*sign(z0_1-f_k)+z2_1)+z1_1;
        z2=h*(-lambda3*L^(3/5)*abs(z0_1-f_k)^(2/5)*sign(z0_1-f_k)+z3_1)+z2_1;
        z3=h*(-lambda4*L^(4/5)*abs(z0_1-f_k)^(1/5)*sign(z0_1-f_k)+z4_1)+z3_1;
        z4=h*(-lambda5*L*sign(z0_1-f_k))+z4_1;
    
	elseif differentiator_order==3
        lambda0=3;lambda1=4.16;lambda2=3.06;lambda3=1.1;        
        z0=h*(-lambda0*L^(1/4)*abs(z0_1-f_k)^(3/4)*sign(z0_1-f_k)+z1_1)+z0_1;    
        z1=h*(-lambda1*L^(2/4)*abs(z0_1-f_k)^(2/4)*sign(z0_1-f_k)+z2_1)+z1_1;
        z2=h*(-lambda2*L^(3/4)*abs(z0_1-f_k)^(1/4)*sign(z0_1-f_k)+z3_1)+z2_1;
        z3=h*(-lambda3*L*sign(z0_1-f_k))+z3_1;
        
    elseif differentiator_order==2
        lambda0=2;lambda1=2.12;lambda2=1.1;        
        z0=h*(-lambda0*L^(1/3)*abs(z0_1-f_k)^(2/3)*sign(z0_1-f_k)+z1_1)+z0_1;    
        z1=h*(-lambda1*L^(2/3)*abs(z0_1-f_k)^(1/3)*sign(z0_1-f_k)+z2_1)+z1_1;
        z2=h*(-lambda2*L*sign(z0_1-f_k))+z2_1;
    
    else
        disp("Error in E-AO-STD")
    end
        

    CT=CT+toc;
    
    if output_order==0
        u=z0;
    elseif output_order==1
        u=z1;
    elseif output_order==2
        u=z2;
    elseif output_order==3
        u=z3;
    elseif output_order==4
        u=z4;     
    elseif output_order==5
        u=z5;     
    elseif output_order==6
        u=z6;     
    elseif output_order==7
        u=z7;          
    else
        disp("Error in E-AO-STD")
    end
    
	calculation_time=CT;    
end

