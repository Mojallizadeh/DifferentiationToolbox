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


function [X1,X2] = high_degree_i_test(lambda1,lambda2,mu,h,f_k,reset_flag,o1,o2)
% Control signal of high-degree STA in explicit way
% Ver 24-11-2019
    acc=0.0001;
    it=10000;
    persistent v_1 x_1

    if reset_flag==1
        x_1 = o1;
        v_1 = o2;
    end

    
	a1=h*lambda1;
    a2=1+2*h^2*lambda2*mu;
    a3=h*lambda1*mu;
    a4=(3/2)*h^2*lambda2*mu^2;
        
    bk=-x_1+f_k-h*v_1;
        
	if  bk<-h^2*lambda2/2
%             f = @(x) a2*x^4 + a4*x^3 + a1*x^2 + a3*x + bk+h^2*lambda2/2;
%             df= @(x)  4*a2*x^3+3*a4*x^2+2*a1*x+a3;
%             [sol,N]= newton(f,df,acc,it,reset_flag);
            [sol]=roots([a4,a3,a2,a1,bk+h^2*lambda2/2]);
            for i=1:length(sol)
               if isreal(sol(i)) && sol(i)>0
                   X11=sol(i);
               end
            end
        X11=X11.^2;
%         pause
	elseif bk>h^2*lambda2/2
%             f = @(x) -a2*x^4 - a4*x^3 - a1*x^2 - a3*x + bk-h^2*lambda2/2;
%             df= @(x)  -4*a2*x^3-3*a4*x^2-2*a1*x-a3;
%             [sol,N]= newton(f,df,acc,it,reset_flag);
            [sol]=roots([-a4,-a3,-a2,-a1,bk-h^2*lambda2/2]);
            for i=1:length(sol)
               if isreal(sol(i)) && sol(i)<0
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
    
    

    
% 	if N>N_max
%        N_max=N; 
% 	end
%     N_total=N_total+N;
    
    X1=x_1;
    X2=u_1;
end

