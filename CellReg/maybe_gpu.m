function x=maybe_gpu(x)
% This function moves an array to the GPU when GPU acceleration is enabled
% (see cellreg_use_gpu) and returns it unchanged otherwise, so that the same
% code can run on either device.

% Inputs:
% 1. x - numeric array

% Outputs:
% 1. x - gpuArray or the original array

if cellreg_use_gpu() && ~isa(x,'gpuArray')
    x=gpuArray(x);
end

end
