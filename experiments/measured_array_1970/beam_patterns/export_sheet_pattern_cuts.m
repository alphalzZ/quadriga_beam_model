function cuts=export_sheet_pattern_cuts(qd_path)
%EXPORT_SHEET_PATTERN_CUTS Sample the actual calibrated QuaDRiGa beam fields.
% Uses sheet angles only; the legacy rectangular cache weights are not used.
if nargin<1,qd_path=getenv('QD_PATH');end
if ~isempty(qd_path),addpath(genpath(qd_path));end
assert(exist('qd_arrayant','file')~=0,'Provide the QuaDRiGa path');
here=fileparts(mfilename('fullpath'));addpath(fullfile(here,'..','time_domain'));
source_cache=fullfile(here,'..','..','..','src','beam_codebook','results','beam_model_sheet.mat');
S=load(source_cache,'ele_tab','azi_tab','beam_name');
[~,array,info]=make_allfit_link_arrays(S.ele_tab,S.azi_tab,'RHCP',0);
theta_deg=unique([0:0.25:90 S.ele_tab]);
phi_deg=unique([0:0.25:360 mod(S.azi_tab,360)]);
K=numel(S.ele_tab);
theta_gain_db=zeros(K,numel(theta_deg));phi_gain_db=zeros(K,numel(phi_deg));
normalization_db=zeros(K,1);peak_error_deg=zeros(K,2);
for k=1:K
    % Matrix coordinates always reference the fixed array normal, not beam axis.
    [v,h]=array.interpolate(ones(size(theta_deg))*S.azi_tab(k)*pi/180, ...
        (90-theta_deg)*pi/180,k);
    theta_gain_db(k,:)=10*log10(max(abs(v).^2+abs(h).^2,realmin));
    [v,h]=array.interpolate(phi_deg*pi/180,ones(size(phi_deg))*(90-S.ele_tab(k))*pi/180,k);
    phi_gain_db(k,:)=10*log10(max(abs(v).^2+abs(h).^2,realmin));
    assert(abs(phi_gain_db(k,1)-phi_gain_db(k,end))<1e-7,'Azimuth seam does not close');
    normalization_db(k)=max([theta_gain_db(k,:) phi_gain_db(k,:)]);
    [~,it]=max(theta_gain_db(k,:));[~,ip]=max(phi_gain_db(k,:));
    peak_error_deg(k,1)=abs(theta_deg(it)-S.ele_tab(k));
    dotp=cosd(S.ele_tab(k))^2+sind(S.ele_tab(k))^2*cosd(phi_deg(ip)-S.azi_tab(k));
    peak_error_deg(k,2)=acosd(max(-1,min(1,dotp)));
end
assert(all(isfinite([theta_gain_db(:);phi_gain_db(:)])));
assert(max(peak_error_deg(:))<2.5,'Cut peaks do not agree with requested beam directions');
cuts=struct('beam_name',{S.beam_name},'pitch_deg',S.ele_tab,'azimuth_deg',S.azi_tab, ...
    'theta_deg',theta_deg,'phi_deg',phi_deg,'theta_gain_db',theta_gain_db, ...
    'phi_gain_db',phi_gain_db,'normalization_db',normalization_db, ...
    'theta_normalized_db',bsxfun(@minus,theta_gain_db,normalization_db), ...
    'phi_normalized_db',bsxfun(@minus,phi_gain_db,normalization_db), ...
    'extrapolated',info.extrapolated,'fc_hz',info.fc,'elements',info.N, ...
    'peak_error_deg',peak_error_deg,'polarization','RHCP', ...
    'angle_convention','fixed-normal theta=sheet ele, phi=sheet azi', ...
    'normalization','one shared maximum of both sampled cuts per beam', ...
    'field_source','same QuaDRiGa calibrated beam interpolator as time-domain link');
out=fullfile(here,'results');if ~exist(out,'dir'),mkdir(out);end
save(fullfile(out,'sheet_pattern_cuts.mat'),'-struct','cuts','-v7');
fprintf('PASS: %d beams; seam closure, finite powers and cut directions; max error %.4f deg\n', ...
    K,max(peak_error_deg(:)));
end
