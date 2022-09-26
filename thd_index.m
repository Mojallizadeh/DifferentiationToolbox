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

function [thd_out] = thd_index(u,h,t_f,w)
% Calculate THD

    fftSignal = fft(u);
    fftSignal = abs(fftshift(fftSignal));
    f = 1/(2*h)*linspace(-1,1,t_f/h+1);


    [~,i] = min(abs(f-w/(2*pi)));

	v0=fftSignal(i);
    
    thd_out=100*(sum(fftSignal)/length(fftSignal))/(v0);

% length(f)
% length(fftSignal)
% figure
%      plot(f,fftSignal)
% %      figure

end

