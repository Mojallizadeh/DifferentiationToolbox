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


function [u,calculation_time] = high_degree_e(lambda1,lambda2,mu,h,f_k,reset_flag,order)
% Control signal of high-degree STA in explicit way
% Ver 24-11-2019  Ok

    persistent v_1 v_2 v_3 x_1 x_2 x_3 CT
    
    if reset_flag==1
        x_1 = 0;
        x_2 = 1;
        x_3 = 0;
        v_1 = 1;
        v_2 = 0;
        v_3 = -1;
        CT=0;
    end
    
    tic
    
        sigma=x_1-f_k;
        
            u_1=v_1-lambda1*sqrt(abs(sigma))*sign(sigma)-lambda1*mu*abs(sigma)^(3/2)*sign(sigma);
            v_1=-h*lambda2*((1/2)*sign(sigma)+2*mu*sigma+(3/2)*mu^2*sigma^2*sign(sigma))+v_1;
            x_1=h*u_1+x_1;
        
        u=u_1;
        
        if order==2 || order==3
        
            sigma=x_2-u_1;
            u_2=v_2-lambda1*sqrt(abs(sigma))*sign(sigma)-lambda1*mu*abs(sigma)^(3/2)*sign(sigma);
            v_2=-h*lambda2*((1/2)*sign(sigma)+2*mu*sigma+(3/2)*mu^2*sigma^2*sign(sigma))+v_2;
            x_2=h*u_2+x_2;
        
        u=u_2;
        
        end
        
        if order==3
        
            sigma=x_3-u_2;
            u_3=v_3-lambda1*sqrt(abs(sigma))*sign(sigma)-lambda1*mu*abs(sigma)^(3/2)*sign(sigma);
            v_3=-h*lambda2*((1/2)*sign(sigma)+2*mu*sigma+(3/2)*mu^2*sigma^2*sign(sigma))+v_3;
            x_3=h*u_3+x_3;
        
            u=u_3;
            
        end
        
        CT=CT+toc;

        calculation_time=CT;
    
end

