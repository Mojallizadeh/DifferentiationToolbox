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


function [u,calculation_time] = kalman_diff(R,h,f,reset_flag,order,x0)
% Control signal of STA differentiato in implicit way
% 24-11-2019 ok.

    persistent P A H Q x CT
    
    if reset_flag==1
        P = zeros(4);
        A=[1 h h^2/2 h^3/6;...
             0  1  h h^2/2;...
             0  0  1 h;...
             0  0  0  1];
         H=[1 0 0 0];

         Q=[0 0 0 0;...
            0 0 0 0;...
            0 0  0 0;...
            0 0 0 h];

        x=[x0(1);x0(2);x0(3);x0(4)];
        
        CT=0;
    end
    tic % Timer

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    x=A * x; 

    P=A*P*A'+Q;

    y=f-H*x;

    S=H*P*H'+R;

    K=P*H'/S;

    x=x+K*y;

    P=P-K*H*P;
        
    
    if order==1
        u=x(2);
    elseif order==2
        u=x(3);
    elseif order==3
        u=x(4);
    else
        u=0.1;
    end
    

	CT=CT+toc;

	calculation_time=CT;

end

