function model = run_beam_model_triangular_cp(opts)
%RUN_BEAM_MODEL_TRIANGULAR_CP Planar lattice with ideal circular-polarized ports.
%   model = run_beam_model_triangular_cp(); % 3.5 GHz, 32x16, triangular, RHCP
%   model = run_beam_model_triangular_cp(struct('polarization','DUAL'));
%   Options: lattice ('triangular'/'rectangular'), polarization ('LHCP',
%   'RHCP','DUAL'), fc, Mv, Mh, spacing_lambda, angle_source ('photo'/'sheet'),
%   grid_step_deg, output_dir. output_dir='' disables exports.
%
% Each spatial site has one CP port, or two co-located independent CP ports.
% Jones vectors in QuaDRiGa's local [Fa;Fb] basis are [1;+j]/sqrt(2) for LHCP
% and [1;-j]/sqrt(2) for RHCP. These ideal fields keep the TOTAL element gain
% at elem_gain_dB, not an extra 3 dB per orthogonal component.
% Angular plots use a matched local incident Jones vector per beam. They
% are receive-response patterns, not a full satellite propagation model.
if nargin < 1, opts = struct(); end
here = fileparts(mfilename('fullpath'));
cfg = struct('fc',3.5e9,'Mv',32,'Mh',16,'spacing_lambda',0.5, ...
    'lattice','triangular','polarization','RHCP','angle_source','photo', ...
    'grid_step_deg',1,'output_dir',here);
assert(isstruct(opts) && isscalar(opts), 'opts must be a scalar struct');
names = fieldnames(opts);
for i = 1:numel(names)
    assert(isfield(cfg,names{i}), ['Unknown option: ' names{i}]);
    cfg.(names{i}) = opts.(names{i});
end
for key = {'fc','spacing_lambda','grid_step_deg','Mv','Mh'}
    v = cfg.(key{1});
    assert(isnumeric(v) && isscalar(v) && isreal(v) && isfinite(v) && v>0, ...
        [key{1} ' must be a positive finite scalar']);
end
assert(cfg.Mv==fix(cfg.Mv) && cfg.Mh==fix(cfg.Mh), 'Mv/Mh must be integers');
assert(cfg.grid_step_deg<=10 && abs(180/cfg.grid_step_deg-round(180/cfg.grid_step_deg))<1e-10, ...
    'grid_step_deg must divide 180 and not exceed 10 degrees');
cfg.lattice = lower(cfg.lattice); cfg.polarization = upper(cfg.polarization);
assert(any(strcmp(cfg.lattice,{'triangular','rectangular'})), 'Invalid lattice');
assert(any(strcmp(cfg.polarization,{'LHCP','RHCP','DUAL'})), 'Invalid polarization');
assert(any(strcmp(cfg.angle_source,{'photo','sheet'})), 'Invalid angle_source');

%% Geometry: rows along z, columns along y, broadside along +x.
fc = cfg.fc; lambda = 299792458/fc; d = lambda*cfg.spacing_lambda;
[col,row] = meshgrid(0:cfg.Mh-1,0:cfg.Mv-1);
if strcmp(cfg.lattice,'triangular')
    y = d*(col+0.5*mod(row,2)); z = d*sqrt(3)/2*row;
else
    y = d*col; z = d*row;
end
y = y-mean(y(:)); z = z-mean(z(:));
p_site = [zeros(1,numel(y)); y(:).'; z(:).'];
Nsite = size(p_site,2);

%% Ideal CP ports. Dual-polarization ordering: all LHCP, then all RHCP.
if strcmp(cfg.polarization,'DUAL')
    pol_labels = {'LHCP','RHCP'};
else
    pol_labels = {cfg.polarization};
end
P = numel(pol_labels); J = zeros(2,P);
for q = 1:P
    sign_j = 1;
    if strcmp(pol_labels{q},'RHCP'), sign_j = -1; end
    J(:,q) = [1;sign_j*1j]/sqrt(2);
end
p = repmat(p_site,1,P); N = size(p,2);  % N denotes signal ports in this cache
port_site = repmat(1:Nsite,1,P);
port_pol = kron(1:P,ones(1,Nsite));
port_jones = J(:,port_pol);

%% Same 14 commanded beams and photo interpretation as run_beam_model.m.
beam_base = arrayfun(@(k) sprintf('C%d',k),1:14,'UniformOutput',false);
ele = [45.43 36.28 27.8 20.39 14.59 11.14 32 41.55 50.37 12.93 2.29 8.82 15.58 23.35];
azi = [188.62 184.66 195.4 191.44 222.23 290.9 357.42 358.16 359.61 185.17 24.74 7.02 6.66 5.33];
if strcmp(cfg.angle_source,'photo')
    px = [1497 1392 1225 1137 884 690 268 178 93 1012 787 688 583 483];
    azi = (px-801)/1602*120;
end
K = numel(ele); B = K*P;
[Ws,Us] = build_beam_codebook(ele,azi,p_site,fc);
W = zeros(N,B); beam_name = cell(1,B);
for q = 1:P
    ports = (q-1)*Nsite+(1:Nsite); beams = (q-1)*K+(1:K);
    W(ports,beams) = Ws;
    for k = 1:K, beam_name{beams(k)} = [beam_base{k} '_' pol_labels{q}]; end
end
U = repmat(Us,1,P); beam_pol = kron(1:P,ones(1,K));
beam_jones = J(:,beam_pol);
ele_tab = repmat(ele,1,P); azi_tab = repmat(azi,1,P);

%% Array response on the same spherical angular grid as the original model.
az_g = -180:cfg.grid_step_deg:180; el_g = -90:cfg.grid_step_deg:90;
nA = numel(az_g); nE = numel(el_g); k0 = 2*pi/lambda;
TG = zeros(K,nE,nA); cut_site = zeros(K,nA);
pk_G = zeros(1,K); pk_az = zeros(1,K); pk_el = zeros(1,K);
bw_az = zeros(1,K); bw_el = zeros(1,K);
for ie = 1:nE
    el = el_g(ie);
    uu = [cosd(el)*cosd(az_g); cosd(el)*sind(az_g); sind(el)*ones(1,nA)];
    af = abs(Ws'*exp(1j*k0*(p_site.'*uu))).^2;
    gain = bsxfun(@plus,10*log10(max(af,realmin)),elem_gain_dB(az_g,el));
    TG(:,ie,:) = reshape(gain,K,1,nA);
end
for k = 1:K
    [~,ie] = min(abs(el_g+ele(k)));
    [~,ia] = min(abs(angle_diff_deg(az_g,azi(k))));
    cut_site(k,:) = reshape(TG(k,ie,:),1,nA);
    surface = reshape(TG(k,:,:),nE,nA);
    [pk_G(k),idx] = max(surface(:)); [i,j] = ind2sub(size(surface),idx);
    pk_az(k) = az_g(j); pk_el(k) = el_g(i);
    % Contiguous 3 dB AF widths around the commanded grid sample, not the
    % potentially mirrored back lobe. Resolution is grid_step_deg.
    af_az = cut_site(k,:)-elem_gain_dB(az_g,el_g(ie));
    af_el = surface(:,ia).'-elem_gain_dB(az_g(ia),el_g);
    bw_az(k) = local_width(af_az,ia,cfg.grid_step_deg);
    bw_el(k) = local_width(af_el,ie,cfg.grid_step_deg);
end
cutAz = repmat(cut_site,P,1);
gain_cmd = 10*log10(Nsite)+elem_gain_dB(azi_tab,-ele_tab);

%% Polarization-aware response for each commanded user direction.
% A(port,j) = spatial phase * <port Jones, incident Jones_j>.
% Opposite ideal CP has zero overlap in this common local transverse basis.
A = exp(1j*k0*(p.'*U)).*(port_jones'*beam_jones);
Gx = abs(W'*A).^2;
ae = elem_gain_dB(azi_tab,-ele_tab);
XT = bsxfun(@minus,bsxfun(@plus,10*log10(Gx/Nsite),ae),ae.');
XTo = XT; XTo(1:B+1:end) = -Inf;
[worst_xt,worst_j] = max(XTo,[],2);
dist_impl = 23.5./tand(ele_tab);
nearest = Inf;
if Nsite > 1
    for a = 1:Nsite
        ds = sum((p_site-p_site(:,a)).^2,1); ds(a) = Inf;
        nearest = min(nearest,sqrt(min(ds)));
    end
end
assert(max(abs(sum(abs(W).^2,1)-1))<1e-12, 'Weights must have unit norm');
assert(max(abs(diag(Gx)-Nsite))<1e-8, 'Commanded AF peak mismatch');
if Nsite>1, assert(abs(nearest-d)<1e-10, 'Nearest-neighbor spacing mismatch'); end
fprintf('%s / %s: %.2f GHz, %dx%d sites, %d ports, %d beams\n', ...
    cfg.lattice,cfg.polarization,fc/1e9,cfg.Mv,cfg.Mh,N,B);
fprintf('Nearest site spacing %.6g m; weight normalization and AF peak verified\n',nearest);

model = struct('cfg',cfg,'fc',fc,'N',N,'Nsite',Nsite,'K',K,'B',B, ...
    'p',p,'p_site',p_site,'port_site',port_site,'port_pol',port_pol, ...
    'port_jones',port_jones,'pol_labels',{pol_labels},'beam_pol',beam_pol, ...
    'beam_jones',beam_jones,'beam_name',{beam_name},'ele_tab',ele_tab, ...
    'azi_tab',azi_tab,'W',W,'U',U,'XT',XT,'Gx',Gx,'az_g',az_g, ...
    'el_g',el_g,'cutAz',cutAz,'gain_cmd',gain_cmd,'nearest_spacing_m',nearest);
if isempty(cfg.output_dir), return; end
if ~exist(cfg.output_dir,'dir'), mkdir(cfg.output_dir); end
prefix = fullfile(cfg.output_dir,['beam_model_' cfg.lattice '_' lower(cfg.polarization)]);
save([prefix '.mat'],'-struct','model','-v7');
fid = fopen([prefix '_summary.csv'],'w'); assert(fid>=0,'Cannot open summary');
closer = onCleanup(@() fclose(fid));
fprintf(fid,'beam,polarization,ele_deg,azi_deg,dist_m,gain_at_cmd_dBi,gmax_dBi,azmax_deg,elmax_deg,bw_az_deg,bw_el_deg,nbr_beam,xtalk_dB\n');
for b = 1:B
    k = mod(b-1,K)+1;
    fprintf(fid,'%s,%s,%.2f,%.2f,%.1f,%.4f,%.4f,%.2f,%.2f,%.2f,%.2f,%s,%.4f\n', ...
        beam_name{b},pol_labels{beam_pol(b)},ele_tab(b),azi_tab(b),dist_impl(b), ...
        gain_cmd(b),pk_G(k),pk_az(k),pk_el(k),bw_az(k),bw_el(k),beam_name{worst_j(b)},worst_xt(b));
end
dlmwrite([prefix '_cuts.csv'],cutAz,'precision','%.6g');
dlmwrite([prefix '_crosstalk.csv'],XT,'precision','%.6g');
dlmwrite([prefix '_sites.csv'],p_site.','precision','%.12g');
fprintf('Saved %s.mat and summary/cuts/crosstalk/sites CSV files\n',prefix);
end

function width = local_width(gain_db,center,step)
threshold = gain_db(center)-3;
l = center; r = center;
while l>1 && gain_db(l-1)>=threshold, l=l-1; end
while r<numel(gain_db) && gain_db(r+1)>=threshold, r=r+1; end
width = (r-l+1)*step;
end
