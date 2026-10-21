function [a, meta] = build_measured_beam1970(pitch, azimuth, peak_gain_dbi, polarization)
%BUILD_MEASURED_BEAM1970 Power-equivalent beam-state surrogate, not an aperture fit.
%   a = build_measured_beam1970(20,30);
% Only one selected steering state and one chosen ideal CP port are represented.
% The other 77 files are other steering states of THE SAME
% antenna, not 77 independent physical elements. Fa/Fb phase is synthesized.
% Optional peak_gain_dbi is an external calibration offset; default 0 is NOT
% a measured absolute gain. Unmeasured rear hemisphere is zero-filled.
if nargin<3, peak_gain_dbi=0; end
if nargin<4, polarization='RHCP'; end
assert(any(strcmpi(polarization,{'LHCP','RHCP'})),'Choose LHCP or RHCP');
sign_j=-1; if strcmpi(polarization,'LHCP'), sign_j=1; end
assert(exist('qd_arrayant','file')~=0,'Add QuaDRiGa quadriga_src to the path');
assert(isscalar(pitch) && any(pitch==0:10:50));
assert(isscalar(azimuth) && any(azimuth==0:30:360));
assert(isscalar(peak_gain_dbi) && isfinite(peak_gain_dbi));
here=fileparts(mfilename('fullpath'));
source=fullfile(here,'..','..','..','data','measured_patterns', ...
    sprintf('pattern_Pitch_%d_Azimuth_%d.txt',pitch,azimuth));
raw=dlmread(source,'\t',6,0);
assert(isequal(size(raw),[91 361]) && all(raw(:,1)==(0:90).'));
gain=raw(:,2:end);
az_deg=-180:180; el_deg=-90:90;
amp=zeros(181,361);
% Preserve the source values exactly at each measured theta/phi grid sample.
for ie=91:181
    theta=90-el_deg(ie);
    amp(ie,:)=10.^((gain(theta+1,mod(az_deg,360)+1)+peak_gain_dbi)/20);
end
a=qd_arrayant([]); a.set_grid(az_deg*pi/180,el_deg*pi/180,0);
a.no_elements=1; a.center_frequency=1.97e9;
a.name=sprintf('measured1970_P%d_A%d_%s',pitch,azimuth,upper(polarization));
a.Fa=amp/sqrt(2);
a.Fb=sign_j*1j*a.Fa;
a.element_position=zeros(3,1); a.coupling=1;
meta=struct('source',source,'pitch',pitch,'azimuth',azimuth, ...
    'calibration_offset_dbi',peak_gain_dbi,'phase_assumption','zero common phase; ideal CP', ...
    'back_hemisphere','unmeasured; zero-filled','measured_polarization',false, ...
    'is_geometry_identification',false,'polarization',upper(polarization));
end
