function data=export_sheet_pattern_grid(qd_path)
%EXPORT_SHEET_PATTERN_GRID Sample all 14 calibrated beams on the source grid.
if nargin<1,qd_path=getenv('QD_PATH');end
if ~isempty(qd_path),addpath(genpath(qd_path));end
assert(exist('qd_arrayant','file')~=0,'Provide the QuaDRiGa path');
here=fileparts(mfilename('fullpath'));addpath(fullfile(here,'..','time_domain'));
source=fullfile(here,'..','..','..','src','beam_codebook','results','beam_model_sheet.mat');
S=load(source,'ele_tab','azi_tab','beam_name');
[~,array,info]=make_allfit_link_arrays(S.ele_tab,S.azi_tab,'RHCP',0);
theta_deg=0:90;phi_deg=0:359;
[theta,phi]=ndgrid(theta_deg,phi_deg);
K=numel(S.ele_tab);power=zeros(91,360,K);peak_error_deg=zeros(1,K);
for k=1:K
    [v,h]=array.interpolate(phi(:).'*pi/180,(90-theta(:).')*pi/180,k);
    power(:,:,k)=reshape(abs(v).^2+abs(h).^2,91,360);
    [peak,idx]=max(reshape(power(:,:,k),[],1));assert(peak>0);
    u=[sind(theta(idx))*cosd(phi(idx));sind(theta(idx))*sind(phi(idx));cosd(theta(idx))];
    peak_error_deg(k)=acosd(max(-1,min(1,dot(u,info.U_canonical(:,k)))));
end
assert(all(isfinite(power(:))) && all(power(:)>=0));
assert(max(peak_error_deg)<2.5);
data=struct('theta_deg',theta_deg,'phi_deg',phi_deg,'power',power, ...
    'beam_name',{S.beam_name},'pitch_deg',S.ele_tab,'azimuth_deg',S.azi_tab, ...
    'extrapolated',info.extrapolated,'peak_error_deg',peak_error_deg, ...
    'fc_hz',info.fc,'elements',info.N,'polarization','RHCP', ...
    'field_source','calibrated QuaDRiGa beam field power, same as time-domain link');
out=fullfile(here,'results');if ~exist(out,'dir'),mkdir(out);end
save(fullfile(out,'sheet_pattern_grid.mat'),'-struct','data','-v7');
fprintf('PASS: 14 complete 91x360 patterns; max peak error %.6f deg\n',max(peak_error_deg));
end
