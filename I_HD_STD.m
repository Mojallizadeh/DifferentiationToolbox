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


function [u,calculation_time] = I_HD_STD(L,mu,h,f_k,reset_flag,order,x0)
% Control signal of high-degree STA in explicit way
% Ver 24-11-2019

    lambda1=sqrt(L)*1.5;
    lambda2=L*1.1;

    persistent v_1 v_2 v_3 x_1 x_2 x_3 CT

    if reset_flag==1
        x_1 = x0(1);
        x_2 = x0(2);
        x_3 = x0(3);
        v_1 = x0(2);
        v_2 = x0(3);
        v_3 = x0(4);
        CT = 0;      
    end
    tic
    
	a1=h*lambda1;
    a2=1+2*h^2*lambda2*mu;
    a3=h*lambda1*mu;
    a4=(3/2)*h^2*lambda2*mu^2;
        
    bk=-x_1+f_k-h*v_1;
        
	if  bk<-h^2*lambda2/2
            [sol]=roots([a4,a3,a2,a1,bk+h^2*lambda2/2]);
            for i=1:length(sol)
               if isreal(sol(i))
                   X11=sol(i);
               end
            end
        X11=X11.^2;
	elseif bk>h^2*lambda2/2
            [sol]=roots([-a4,-a3,-a2,-a1,bk-h^2*lambda2/2]);
            for i=1:length(sol)
               if isreal(sol(i))
                   X11=sol(i);
               end
            end
            
            X11=-X11.^2;
%             pause
    else
        N=0;
           
	end
    
    if bk<-h^2*lambda2/2
        v_1=-(h*lambda2/2)-2*h*lambda2*mu*X11-(3/2)*h*lambda2*mu^2*X11^2+v_1;
        u_1=v_1-lambda1*sqrt(abs(X11))-lambda1*mu*abs(X11)^(3/2);
    elseif bk>h^2*lambda2/2
        v_1=+(h*lambda2/2)-2*h*lambda2*mu*X11+(3/2)*h*lambda2*mu^2*X11^2+v_1;
        u_1=v_1+lambda1*sqrt(abs(X11))+lambda1*mu*abs(X11)^(3/2);
    else
        v_1=-(x_1-f_k)/h;
        u_1=v_1;
        
    end

    x_1=x_1+h*u_1;

    u=u_1;
    
    
    if order==2 || order==3
       

    bk=-x_2+u_1-h*v_2;
        
	if  bk<-h^2*lambda2/2
            [sol]=roots([a2,a4,a1,a3,bk+h^2*lambda2/2]);
            for i=1:length(sol)
               if isreal(sol(i))
                   X11=sol(i);
               end
            end
        X11=X11.^2;
	elseif bk>h^2*lambda2/2
            [sol]=roots([-a2,-a4,-a1,-a3,bk-h^2*lambda2/2]);
            for i=1:length(sol)
               if isreal(sol(i))
                   X11=sol(i);
               end
            end
            
            X11=-X11.^2;
           
	end
    
    if bk<-h^2*lambda2/2
        v_2=-(h*lambda2/2)-2*h*lambda2*mu*X11-(3/2)*h*lambda2*mu^2*X11^2+v_2;
        u_2=v_2-lambda1*sqrt(abs(X11))-lambda1*mu*abs(X11)^(3/2);
    elseif bk>h^2*lambda2/2
        v_2=+(h*lambda2/2)-2*h*lambda2*mu*X11+(3/2)*h*lambda2*mu^2*X11^2+v_2;
        u_2=v_2+lambda1*sqrt(abs(X11))+lambda1*mu*abs(X11)^(3/2);
    else
        v_2=-(x_2-u_1)/h;
        u_2=v_2;
        
    end

    x_2=x_2+h*u_2;
    
    u=u_2;
    
    end
    if order==3
        
    bk=-x_3+u_2-h*v_3;
        
	if  bk<-h^2*lambda2/2
            [sol]=roots([a2,a4,a1,a3,bk+h^2*lambda2/2]);
            for i=1:length(sol)
               if isreal(sol(i))
                   X11=sol(i);
               end
            end
        X11=X11.^2;
	elseif bk>h^2*lambda2/2
            [sol]=roots([-a2,-a4,-a1,-a3,bk-h^2*lambda2/2]);
            for i=1:length(sol)
               if isreal(sol(i))
                   X11=sol(i);
               end
            end
            
            X11=-X11.^2;
           
	end
    
    if bk<-h^2*lambda2/2
        v_3=-(h*lambda2/2)-2*h*lambda2*mu*X11-(3/2)*h*lambda2*mu^2*X11^2+v_3;
        u_3=v_3-lambda1*sqrt(abs(X11))-lambda1*mu*abs(X11)^(3/2);
    elseif bk>h^2*lambda2/2
        v_3=+(h*lambda2/2)-2*h*lambda2*mu*X11+(3/2)*h*lambda2*mu^2*X11^2+v_3;
        u_3=v_3+lambda1*sqrt(abs(X11))+lambda1*mu*abs(X11)^(3/2);
    else
        v_3=-(x_3-u_2)/h;
        u_3=v_3;
        
    end

    x_3=x_3+h*u_3;
    
    u=u_3;
 
    end

	CT=CT+toc;
	calculation_time=CT;
    
end

