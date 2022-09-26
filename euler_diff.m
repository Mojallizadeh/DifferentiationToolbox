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



function [u,calculation_time] = euler_diff(h,f_k,reset_flag,order)
% Euler differentiator

    persistent f_k_minus_1 df_minus_1 ddf_minus_1 CT
    
	if reset_flag==1
        CT = 0;
        f_k_minus_1 = 0;
        df_minus_1 = 1;
        ddf_minus_1 = 0;
	end
    tic
    
	df=(f_k-f_k_minus_1)/h;
    u=df;
	f_k_minus_1=f_k;
    
    if order==2 || order==3
        ddf=(df-df_minus_1)/h;
        df_minus_1=df;
        u=ddf;
    end
    
	if order==3
        dddf=(ddf-ddf_minus_1)/h;
        u=dddf;
        ddf_minus_1=ddf;
    end
    
	CT=CT+toc;
    
	calculation_time=CT; 

end

