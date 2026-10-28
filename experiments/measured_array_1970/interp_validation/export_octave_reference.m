function export_octave_reference(quadriga_src)
%EXPORT_OCTAVE_REFERENCE Save make_allfit_link_arrays output as a parity reference.
% Usage: octave --eval "export_octave_reference('/path/to/quadriga/quadriga_src')"
here=fileparts(mfilename('fullpath'));
addpath(genpath(quadriga_src));
addpath(fullfile(here,'..','time_domain'));
[~,~,info]=make_allfit_link_arrays(45.43,188.62,'RHCP');
out=fullfile(here,'results','octave_reference_C1.mat');
power_grid=info.power_grid; az_grid_deg=info.az_grid_deg; el_grid_deg=info.el_grid_deg;
commands=info.corner_commands{1}; weights=info.corner_weights{1};
save('-v7',out,'power_grid','az_grid_deg','el_grid_deg','commands','weights');
printf('saved %s\n',out);
end
