function a = make_cp_arrayant(model,step_deg)
%MAKE_CP_ARRAYANT Build a QuaDRiGa array from a triangular/rectangular CP cache.
%   a = make_cp_arrayant(model); % 5-degree FIELD-pattern sampling by default
%   a = make_cp_arrayant(model,1); % finer sampling, more memory
% Requires QuaDRiGa on the path. Geometry and port order are copied exactly.
% Fa/Fb are synthesized ideal CP fields, not measured crossed-dipole patterns.
% Each co-located CP port is stored as one effective element, identity coupling.
if nargin<2, step_deg=5; end
assert(isnumeric(step_deg) && isscalar(step_deg) && isfinite(step_deg) && ...
    step_deg>0 && step_deg<=10 && abs(90/step_deg-round(90/step_deg))<1e-10, ...
    'step_deg must divide 90 and be between 0 and 10');
assert(exist('qd_arrayant','file')~=0,'Add QuaDRiGa quadriga_src to the path first');
az = -180:step_deg:180; el = -90:step_deg:90;
[AZ,EL] = meshgrid(az,el);
amp = sqrt(10.^(elem_gain_dB(AZ,EL)/10));
a = qd_arrayant([]);
a.set_grid(az*pi/180,el*pi/180,0);
a.no_elements = model.N;
a.center_frequency = model.fc;
a.name = [model.cfg.lattice '_' model.cfg.polarization];
a.Fa = bsxfun(@times,amp,reshape(model.port_jones(1,:),1,1,[]));
a.Fb = bsxfun(@times,amp,reshape(model.port_jones(2,:),1,1,[]));
a.element_position = model.p;
a.coupling = eye(model.N);
end
