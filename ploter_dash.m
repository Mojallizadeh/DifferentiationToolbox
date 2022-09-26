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

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%

function [] = ploter_dash(x,y,plot_type,color)
% Select the plot type based on the parameters

    if plot_type == "N"
        plot(x,y,'color',color,'LineWidth',2,'LineStyle','--')
    elseif plot_type == "Y"
        semilogy(x,abs(y),'color',color,'LineWidth',2,'LineStyle','--')
    elseif plot_type == "X"
        semilogx(x,y,'color',color,'LineWidth',2,'LineStyle','--')
    elseif plot_type == "XY"
        loglog(x,abs(y),'color',color,'LineWidth',2,'LineStyle','--')
    end

end

