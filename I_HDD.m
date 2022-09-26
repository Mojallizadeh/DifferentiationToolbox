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


function [u,calculation_time] = I_HDD(f_k,h,L,reset_flag,output_order,x0)
% Arbitrary differentiator order  
    persistent z0 z1 z2 z3 CT

    if reset_flag==1
        z0=x0(1);z1=x0(2);z2=x0(3);z3=x0(4);CT=0;
    end
    tic
    
        lambda0=3;lambda1=4.16;lambda2=3.06;lambda3=1.1;        
        a1=h*lambda0*L^(1/4);
        a2=h^2*lambda1*L^(2/4);
        a3=h^3*lambda2*L^(3/4)    *   3/2;
        a4=h^4*lambda3*L          *   13/6;
    
        bk=-z0-h*z1-3*h^2*z2/2-13*h^3*z3/6+f_k;
    
        if bk<-a4
           sol=roots([1  a1  a2   a3    a4+bk]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=abs(XX);

        elseif bk>a4
           sol=roots([-1  -a1  -a2   -a3    bk-a4]);
            for i=1:length(sol)
                if isreal(sol(i))
                    XX=sol(i);
                end
            end
            X=-abs(XX);

        else
            X=0;
        end

        s=X;

        if bk<-a4

            z3=h*(-lambda3*L)+z3;
            z2=h*(-lambda2*L^(3/4)*abs(s)^1 +z3  )+z2 ;    
            z1=h*(-lambda1*L^(2/4)*abs(s)^2 +z2  )+z1   +h^2*z3/2;     
            z0=h*(-lambda0*L^(1/4)*abs(s)^3 +z1   )+z0  +h^2*z2/2 +   h^3*z3/6;
        elseif bk>a4
            z3=h*(lambda3*L)+z3;
            z2=h*(lambda2*L^(3/4)*abs(s)^1 +z3  )+z2;        
            z1=h*(lambda1*L^(2/4)*abs(s)^2 +z2  )+z1  +h^2*z3/2;        
            z0=h*(lambda0*L^(1/4)*abs(s)^3 +z1  )+z0  +h^2*z2/2 +   h^3*z3/6;     
        else
            z3=z3+(bk)/h^3;
            z2=z2+h*z3;
            z1=z1+h*z2  +  h^2*z3/2; 
            z0=z0+h*z1   +h^2*z2/2 +   h^3*z3/6;
        end
        
   
   

    if output_order==0
        u=z0;
    elseif output_order==1
        u=z1;
    elseif output_order==2
        u=z2;
    elseif output_order==3
        u=z3; 
    else
        disp("Error in I-HDD")
    end

    CT=CT+toc;
    
	calculation_time=CT;    
end

