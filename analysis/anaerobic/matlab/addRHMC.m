function addRHMC()
% addRHMC  Put PolytopeSamplerMatlab on the path (its initSampler adds the
%   code and compiled-binary folders).
root = getenv('POLYTOPESAMPLER_ROOT');   % checkout of PolytopeSamplerMatlab
assert(isfolder(root), 'Set POLYTOPESAMPLER_ROOT to the PolytopeSamplerMatlab checkout');
here = pwd; cleanup = onCleanup(@() cd(here));
cd(root);
initSampler;
end
