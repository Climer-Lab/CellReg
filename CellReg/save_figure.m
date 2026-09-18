function save_figure(fig,figures_directory,figure_name)
% This function saves a CellReg figure as both .fig and .png without
% crashing the pipeline when the figure can no longer be found (e.g. the
% user closed the window while a later stage was still running, or the
% figures directory disappeared). In those cases a warning is issued and
% the registration continues.

% Inputs:
% 1. fig - figure handle
% 2. figures_directory
% 3. figure_name - file name without extension

if nargin<3 || isempty(figures_directory)
    warning('CellReg:figureNotSaved','No figures directory given - "%s" was not saved.',figure_name)
    return
end

if ~isgraphics(fig,'figure')
    warning('CellReg:figureNotFound','Figure "%s" could not be found (was the window closed?) - it was not saved.',figure_name)
    return
end

try
    set(fig,'PaperPositionMode','auto')
    savefig(fig,fullfile(figures_directory,[figure_name '.fig']))
    saveas(fig,fullfile(figures_directory,figure_name),'png')
catch err
    warning('CellReg:figureNotSaved','Figure "%s" could not be saved (%s) - continuing.',figure_name,err.message)
end

end
