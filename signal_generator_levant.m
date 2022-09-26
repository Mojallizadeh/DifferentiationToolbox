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

function [f,df,ddf,dddf] = signal_generator_levant(len,w,amp,phi,SNR,h,type,An,Wn)
%% Generate the input signal and its first three differentiation
% A sinusoidal signal is considered to investigate the THD effectively
    f_t=zeros(1,len+1);
    f=zeros(1,len+1);
    df=zeros(1,len+1);
    ddf=zeros(1,len+1);
    dddf=zeros(1,len+1);
    
    for k=1:floor(len+1)
        
%         f_t(k)=(h*k)^4-5*(h*k)^2+2*(h*k); % Signal
%         df(k)=4*(h*k)^3-10*(h*k)+2; % Real 1st-order differentiation
%         ddf(k)=12*(h*k)^2-10; % Real 2nd-order differentiation
%         dddf(k)=24*(h*k); % Real 3rd-order differentiation

            f_t(k)=(h*k)+2;
            df(k)=1;
            ddf(k)=0;
            dddf(k)=0;

    end
    

        f=f_t;

    
end

