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

function [r1,r2,im1,im2]=solveq(u,v,n,a);
% Solves x^2 + ux + v = 0 (n ¹ 1) or x + a(1) = 0 (n = 1).
%
% Example call: [r1,r2,im1,im2]=solveq(u,v,n,a)
% r1, r2 are real parts of the roots,
% im1, im2 are the imaginary parts of the roots.
% Called by function bairstow.
%
if n==1
  r1=-a(1);im1=0; r2=0; im2=0;
else
  d=u*u-4*v;
  if d<0
    d=-d;
    im1=sqrt(d)/2; r1=-u/2; r2=r1; im2=-im1;
  elseif d>0
    r1=(-u+sqrt(d))/2; im1=0; r2=(-u-sqrt(d))/2; im2=0;
  else
    r1=-u/2; im1=0; r2=-u/2; im2=0;
  end;
end;