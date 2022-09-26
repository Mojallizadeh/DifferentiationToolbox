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


function [X,N]= halley(f,f_inv_1,f_inv_2,tol,it,reset_flag)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here 


persistent x


    if reset_flag==1
        x=0;
        
        
    end

    e=inf;
    n=0;
    x=0;
    while abs(e)>tol && n<it
        x=x+2* double(subs(f_inv_1,x)) /double(subs(f_inv_2,x));
        e=f(x);
        n=n+1;
    end

    X=x;
    N=n;
end

