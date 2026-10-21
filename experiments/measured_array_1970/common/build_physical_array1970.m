function [a, meta] = build_physical_array1970(beam_index, model_path, step_deg, polarization)
%BUILD_PHYSICAL_ARRAY1970 Fitted triangular aperture, ~500 elements per RX/TX.
% Add QuaDRiGa quadriga_src to the path first.
%   a = build_physical_array1970(1); % one beamformed RHCP signal port
%   a = build_physical_array1970([]); % ~500 independent RHCP element ports
%   a = build_physical_array1970(1,[],5,'LHCP'); % choose opposite CP state
% Build TX and RX arrays separately. Polarization does NOT double either count.
% The array normal is +z. Observation theta=90-elevation, phi=azimuth.
% Element gain scale is ASSUMED 0 dBi; normalized data cannot determine it.
if nargin<1, beam_index=1; end
if nargin<2 || isempty(model_path)
    model_path=fullfile(fileparts(mfilename('fullpath')),'..','baseline_fit','results','ellipse_physical_model.mat');
end
if nargin<3, step_deg=1; end
if nargin<4, polarization='RHCP'; end
assert(any(strcmpi(polarization,{'LHCP','RHCP'})),'Choose LHCP or RHCP');
sign_j=-1; if strcmpi(polarization,'LHCP'), sign_j=1; end
assert(exist('qd_arrayant','file')~=0,'Add QuaDRiGa quadriga_src to the path');
assert(isscalar(step_deg) && step_deg>0 && step_deg<=10 && mod(90,step_deg)==0);
M=load(model_path); N=M.Nsite;
assert(M.fc==1.97e9 && size(M.p_site,2)==N && numel(M.taper)==N);
% The element has no azimuth dependence; three azimuth nodes suffice. Array
% phase is evaluated from positions, not interpolated on this coarse grid.
az=[-pi 0 pi]; el=(-90:step_deg:90)*pi/180;
% Cosine^q power pattern on the measured front hemisphere; no backside data.
amp=repmat(max(sin(el(:)),0).^(M.element_cos_power/2),1,numel(az));
a=qd_arrayant([]); a.set_grid(az,el,0); a.no_elements=N;
a.name=['1970MHz_fitted_triangular_' upper(polarization)]; a.center_frequency=M.fc;
a.element_position=M.p_site;
a.Fa=repmat(amp/sqrt(2),1,1,N);
a.Fb=sign_j*1j*a.Fa;
if isempty(beam_index)
    a.coupling=eye(N);
else
    assert(isscalar(beam_index) && beam_index==fix(beam_index) && ...
        beam_index>=1 && beam_index<=size(M.commands,1),'Invalid command index');
    % QuaDRiGa TX field uses exp(+j*k*p.u); steering therefore uses conjugate
    % of receive weights. Do not conflate receive W^H with TX coupling W.
    w=conj(M.W_receive(:,beam_index));
    a.coupling=w;
end
meta=struct('fc',M.fc,'Nsite',N,'Nports',N,'commands',M.commands, ...
    'beam_index',beam_index,'polarization',upper(polarization), ...
    'element_gain_assumption_dbi',0,'measured_complex_phase',false, ...
    'measured_polarization',false,'measured_absolute_gain',false);
end
