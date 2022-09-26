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


function [u,calculation_time] = alien(T,k,mu,h,f,reset_flag,order,x0)
% Implementation of alien differentiator
% CT=calculation time
    persistent y j CT
    n=order;
    if reset_flag==1
        j=0;
        CT=0;
    end
    tic

    x=0;
    j=j+1;
    
	y(j)=f;
    
    g=factorial(k+mu+2*n+1)/(factorial(k+n)*factorial(mu+n));
     
	if j<=floor(T/h)
        m=j-1;
	else
        m=floor(T/h);
    end
        
    if order==1    
    
        for i=1:m
            tau=i/m;
            xx=tau^(k + n - 1)*(k + n)*(1 - tau)^(mu + n) - tau^(k + n)*(mu + n)*(1 - tau)^(mu + n - 1);
            x=x+(g/T)*(h/T)*(xx)*y(j-i);
        end
    
    elseif order==2
            
        for i=1:m
            tau=i/m;
            xx=tau^(k + n - 2)*(k + n)*(1 - tau)^(mu + n)*(k + n - 1) - 2*tau^(k + n - 1)*(k + n)*(mu + n)*(1 - tau)^(mu + n - 1) + tau^(k + n)*(mu + n)*(1 - tau)^(mu + n - 2)*(mu + n - 1);
            x=x+(g/T)*(h/T^2)*(xx)*y(j-i);
        end
            
    elseif order==3
        
        for i=1:m
            tau=i/m;
            xx=tau^(k + n - 3)*(k + n)*(1 - tau)^(mu + n)*(k + n - 1)*(k + n - 2) - 3*tau^(k + n - 2)*(k + n)*(mu + n)*(1 - tau)^(mu + n - 1)*(k + n - 1) + 3*tau^(k + n - 1)*(k + n)*(mu + n)*(1 - tau)^(mu + n - 2)*(mu + n - 1) - tau^(k + n)*(mu + n)*(1 - tau)^(mu + n - 3)*(mu + n - 1)*(mu + n - 2);
            x=x+(g/T)*(h/T^2)*(xx)*y(j-i);
        end
        
    else
       disp('Error! ALIEN is not defined for n>3'); 
    end
    
    
    
    CT=CT+toc;

    u=x;

    calculation_time=CT;
end

