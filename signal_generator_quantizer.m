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

function [f,df,ddf,dddf] = signal_generator_quantizer(len,w,amp,phi,SNR,h,type,An,Wn)
%% Generate the input signal and its first three differentiation
% A sinusoidal signal is considered to investigate the THD effectively
    f_t=zeros(1,len+1);
    f=zeros(1,len+1);
    df=zeros(1,len+1);
    ddf=zeros(1,len+1);
    dddf=zeros(1,len+1);
    
    for k=1:floor(len+1)
        
        f_t(k)=quant(amp*sin(w*h*k+phi),0.1); % Signal
        df(k)=amp*w*cos(w*h*k+phi); % Real 1st-order differentiation
        ddf(k)=-amp*w^2*sin(w*h*k+phi); % Real 2nd-order differentiation
        dddf(k)=-amp*w^3*cos(w*h*k+phi); % Real 3rd-order differentiation
    end
    
    
    %% Add noise
    if type=="white"
        f=awgn(f_t,SNR,'measured');
    elseif type=="sin"
        for k=1:floor(len+1)
            f(k)=f_t(k)+An*sin(Wn*h*k);
        end
    else
        disp('Undefined noise')
    end
    
end

