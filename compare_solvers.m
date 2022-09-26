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


clc;clear;close all;


P=[1 1 1 1 1 -1];
syms x
f = @(x) x^4+x^3+x^2+x-1;

f_inv=@(x) 1/(x^4+x^3+x^2+x-1);

f_inv_1=diff(f_inv,x);

f_inv_2=diff(f_inv_1,x);

df=@(x) 4*x^3+3*x^2+2*x+1;
tol=1e-5;
it=1e6;
reset_flag=1;

tic


X=roots(P);


MATLAB_real_time=toc


tic
[X,N]= newton(f,df,tol,it,reset_flag);
Newton_real_time=toc

tic
[X,N]= halley(f,f_inv_1,f_inv_2,tol,it,reset_flag);
Halley_real_time=toc

tic
[rts,it]=bairstow(P(2:end),5,tol);
Halley_real_time=toc







