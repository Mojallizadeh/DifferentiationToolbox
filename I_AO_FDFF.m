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


function [u,calculation_time] = I_AO_FDFF(F,eps,ws,wf,a1,ro,f_k,h,reset_flag,order,x0)
%   Quadratic differentiator
%   Detailed explanation goes here

    persistent z1 z2 z3 z4 w2 w3 w4 CT B C
    
    y=f_k;
    
    k1=2.613*wf;
    k2=3.414*wf^2;
    k3=2.613*wf^3;
    k4=1*wf^4;
    
    a2=2*a1*ws;
    a3=2*a1*ws^2;
    a4=a1*ws^3;

    if reset_flag==1
        z1=x0(1);
        z2=x0(2);
        z3=x0(3);
        z4=x0(4);
        w2=x0(2);
        w3=x0(3);
        w4=x0(4);
        C=1+h*k1+h^2*k2+h^3*k3+h^4*k4;
        B=(h*a1+h^2*a2+h^3*a3+h^4*a4)/C;
        CT= 0;
    end
    tic

    z1m=z1;    
    z2m=z2; 
    z3m=z3; 
    z4m=z4;
    
    A=(z1m+h*(z2m+y*k1+h*(z3m...
        +y*k2)+h^2*(z4m+y*k3)...
        +h^3*y*k4))/C;

    z1=B* psi_kikk((y-A)/B,ro)+A;
    
    z2=(a4*h^2*z1-a4*h^2*z1m+a1*z2m...
        +a1*h*z3m+a1*h^2*z4m-a4*h^3*y*k1...
        +a4*h^3*z1*k1+a2*(z1-z1m-h*y*k1...
        +h*z1*k1)+a3*h*(z1-z1m-h*y*k1...
        +h*z1*k1)+a1*h*y*k2-a1*h*z1*k2...
        +a1*h^2*y*k3-a1*h^2*z1*k3+a1*h^3*y*k4...
        -a1*h^3*z1*k4)...
        /(a1+h*(a2+h*(a3+a4*h)));
    
    z3=(a3*(-z1m+z1*(1+h*k1+h^2*k2)...
        -h*(z2m+y*k1+h*y*k2))...
        +a4*h*(-z1m+z1*(1+h*k1+h^2*k2)...
        -h*(z2m+y*k1+h*y*k2))+(a1...
        +a2*h)*(z3m+h*(z4m+(y...
        -z1)*(k3+h*k4))))...
        /(a1+h*(a2+h*(a3+a4*h)));
    
    z4=(-a4*(z1m-z1*(1+h*k1+h^2*k2+h^3*k3)...
        +h*(z2m+h*z3m+y*k1+h*y*k2...
        +h^2*y*k3))+(a1+h*(a2+h*a3))*(z4m...
        +h*(y-z1)*k4))...
        /(a1+h*(a2+h*(a3+a4*h)));
        
    w2=h*F*sat((h*eps*z2+w2)/(h*F*(h*eps+1))-(w2)/(h*F))+w2;
    w3=h*F*sat((h*eps*z3+w3)/(h*F*(h*eps+1))-(w3)/(h*F))+w3;
    w4=h*F*sat((h*eps*z4+w4)/(h*F*(h*eps+1))-(w4)/(h*F))+w4;
    

    out1=w2;
    out2=w3;
    out3=w4;
    
    if order==1
        u=out1;  
    elseif order==2
        u=out2; 
    elseif order==3
        u=out3;          
    else
        u=0;
    end
	    
	CT=CT+toc;
    
	calculation_time=CT;    

end

