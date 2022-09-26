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

function [u,calculation_time,mag,gamma1,alpha1] = VGED(mu,tau,wc,q,h,f_k,reset_flag,order,x0)
% VGED
% Ver 01-06-2020  Ok


    persistent CT j z1 z2 alpha gamma num den f
    
    
    k1=1.5;
    k2=1.1;
    eps=1/mu;
    
    if reset_flag==1
        z1=x0(1);
        z2=x0(2);
        gamma=0;
        CT=0;    
        alpha=0.55;
        j=0;
        s=tf('s');
        s=s/wc;
        G=(s^4)/((s^2+0.7654*s+1)*(s^2+1.8478*s+1));
        Gd=c2d(G,h);
        [num,den]=tfdata(Gd);
        f=0;
    end
    
	if order==1
    
    j=j+1;

    tic
    
    f(j)=f_k;

    out=filter(num{:},den{:},f);
    
    mag=abs(out(j));
    
    
    sigma=z1-f(j);
    
    if sigma==0
       sigma=1e-12; 
    end
    
    
    
    
    z1b=h*(z2-k1*mu*abs(sigma)^(alpha)*sign(sigma))+z1;
    z2b=h*(-k2*alpha*mu^2*(abs(sigma))^(2*alpha-1)*sign(sigma))+z2;
    gamma=h*(-tau*gamma+tau*mag)+gamma;
    alpha=(1/2)*(1+(gamma^q)/(gamma^q+eps));
    
    if alpha>=1
       alpha=0.99; 
    end
    
    
    z1=z1b;
    z2=z2b;
    

        u=z2;
    else
        u=0.1;
        disp("VGED is not designed for higher-orders");
    end
     
    CT=CT+toc;
    gamma1=gamma;
    alpha1=alpha;

    calculation_time=CT;
end
