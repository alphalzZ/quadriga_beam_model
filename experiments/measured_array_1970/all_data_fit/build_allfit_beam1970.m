function [a,meta]=build_allfit_beam1970(pitch,azimuth,polarization,gain_offset_dbi)
%BUILD_ALLFIT_BEAM1970 All-data calibrated beam output for one scan state.
%   a=build_allfit_beam1970(20,30); % default RHCP
%   a=build_allfit_beam1970(20,30,'LHCP');
% All 78 states and every measured theta/phi value enter the residual fit.
% This is ONE equivalent combined-beam signal port. The baseline aperture has
% 502 physical RX elements; a.no_elements=1 does not describe its hardware.
% Residuals are state/direction dependent, not recovered element-level weights.
% Phase comes from the assumed physical model; polarization and absolute gain
% are also assumptions. Unmeasured rear hemisphere is zero-filled.
if nargin<3,polarization='RHCP';end
if nargin<4,gain_offset_dbi=0;end
assert(exist('qd_arrayant','file')~=0,'Add QuaDRiGa quadriga_src to the path');
assert(any(strcmpi(polarization,{'LHCP','RHCP'})),'Choose LHCP or RHCP');
assert(isscalar(pitch) && isscalar(azimuth) && isscalar(gain_offset_dbi) && isfinite(gain_offset_dbi));
path=fullfile(fileparts(mfilename('fullpath')),'results','calibration_tables.mat');
info=dir(path);assert(~isempty(info),'Run fit_all_data.py first');
persistent tables cached_stamp
stamp=[info.datenum info.bytes];
if isempty(tables) || ~isequal(cached_stamp,stamp)
    tables=load(path);cached_stamp=stamp;
end
index=find(tables.commands(:,1)==pitch & tables.commands(:,2)==azimuth);
assert(numel(index)==1,'Select a measured Pitch (0:10:50) / Azimuth (0:30:360)');
gain=tables.baseline_db(:,:,index)+tables.residual_db(:,:,index);
phase=double(tables.assumed_phase_rad(:,:,index));
field=10.^((gain+gain_offset_dbi)/20).*exp(1j*phase);
az_deg=-180:180;el_deg=-90:90;
effective=zeros(181,361);
for ie=91:181
    theta=90-el_deg(ie);
    effective(ie,:)=field(theta+1,mod(az_deg,360)+1);
end
a=qd_arrayant([]);a.set_grid(az_deg*pi/180,el_deg*pi/180,0);
a.center_frequency=1.97e9;a.no_elements=1;
a.name=sprintf('allDataFit1970_P%d_A%d_%s',pitch,azimuth,upper(polarization));
a.Fa=effective/sqrt(2);
sign_j=-1;if strcmpi(polarization,'LHCP'),sign_j=1;end
a.Fb=sign_j*1j*a.Fa;
a.element_position=zeros(3,1);a.coupling=1;
meta=struct('pitch',pitch,'azimuth',azimuth,'fit_files',78,'holdout_files',0, ...
    'physical_receive_elements',502,'equivalent_signal_ports',1, ...
    'is_element_level_model',false,'polarization',upper(polarization), ...
    'phase_source','assumed physical baseline, not measurement', ...
    'measured_absolute_gain',false,'measured_polarization',false, ...
    'gain_offset_dbi',gain_offset_dbi,'calibration_cache',path);
end
