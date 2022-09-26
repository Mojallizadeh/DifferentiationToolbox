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

function [X1,X2] = STA_implicit_test(lambda1,lambda2,h,f_k,reset_flag,o1,o2)
% Control signal of STA differentiato in implicit way
% 24-11-2019 ok.

    persistent v_1 x_1 u_1
    
    if reset_flag==1
        x_1 =o1;
        v_1 =o2;
    end
    
    a=h*lambda1; % a
    b_k=-x_1-h*v_1+f_k; % b
    
	% Case 1: ----------------------------------------------------
	if b_k<-h^2*lambda2    
        sqrt_abs_x_tilde_next=(-a+sqrt(a^2-4*(b_k+lambda2*h^2)))/(2);
        v_1=v_1-h*lambda2;
        u_1=-lambda1*sqrt_abs_x_tilde_next+v_1;
            

    % Case 2: ----------------------------------------------------   
	elseif (b_k>-h^2*lambda2  && b_k<h^2*lambda2)
        v_1=v_1+b_k/h;
        u_1=-(x_1-f_k)/h;

    % Case 3: ---------------------------------------------------- 
    else
        sqrt_abs_x_tilde_next=(-a+sqrt(a^2+4*(b_k-lambda2*h^2)))/(2);
        v_1=v_1+h*lambda2;
        u_1=lambda1*sqrt_abs_x_tilde_next+v_1;

	end
        
    x_1=h*(u_1)+x_1;
    
    X1=x_1;
    X2=u_1;
    
    
	end
    

