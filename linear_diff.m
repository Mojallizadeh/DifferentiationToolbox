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


function [out,calculation_time] = LD(f,c,h,reset_flag,order)
% Linear differentiator

    persistent y1 y2 y3 f_11 f_12 f_13 CT
    
    if reset_flag==1
        y1 =1;
        y2=0;
        y3=-1;
        f_11=0;
        f_12=1;
        f_13=0;        
        CT=0;
    end
    
    tic % Start the timer
    
    y1=(y1+c*(f-f_11))/(1+h*c);
	f_11=f;
    out1=y1;
    
    if order==2 || order==3
        y2=(y2+c*(y1-f_12))/(1+h*c);
        f_12=y1;
        out1=y2;    
    end
    
    if order==3
        y3=(y3+c*(y2-f_13))/(1+h*c);
        f_13=y2;
        out1=y3;    
    end
    
	CT=CT+toc;
	calculation_time=CT;  
    out=out1;

end
