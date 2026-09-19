function use_gpu=cellreg_use_gpu(setting)
% This function stores and returns whether GPU acceleration should be used
% for the heavy image operations (currently the rotation search and the
% cross-correlations in align_images, and the non-rigid registration).
%
%   cellreg_use_gpu()        - returns the current setting. On the first call
%                              it is initialized to true if a supported GPU
%                              and the Parallel Computing Toolbox are
%                              available, false otherwise.
%   cellreg_use_gpu(true)    - request GPU use (falls back to false with a
%                              warning if no usable GPU is found)
%   cellreg_use_gpu(false)   - disable GPU use
%   cellreg_use_gpu('auto')  - re-detect

% Inputs:
% 1. setting (optional) - true / false / 'auto'

% Outputs:
% 1. use_gpu - logical

persistent current_setting

if nargin>0
    if ischar(setting) || isstring(setting)
        current_setting=detect_gpu();
    elseif setting
        if detect_gpu()
            current_setting=true;
        else
            warning('CellReg: no usable GPU was found - continuing on the CPU');
            current_setting=false;
        end
    else
        current_setting=false;
    end
elseif isempty(current_setting)
    current_setting=detect_gpu();
end
use_gpu=current_setting;

end

function available=detect_gpu()
available=false;
try
    if license('test','Distrib_Computing_Toolbox') && exist('canUseGPU','file') && canUseGPU()
        available=true;
    end
catch
    available=false;
end
end
