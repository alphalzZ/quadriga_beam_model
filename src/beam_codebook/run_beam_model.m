% RUN_BEAM_MODEL  Model the 14 receive beams given in the spreadsheet (GNU Octave)
%   Requires build_beam_codebook.m, elem_gain_dB.m, angle_diff_deg.m in the same folder.
more off;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(here);
addpath(fullfile(here,'..','common'));
out=fullfile(here,'results');if ~exist(out,'dir'),mkdir(out);end

%% ---- inputs -------------------------------------------------------------
fc      = 3.5e9;          % carrier frequency [Hz]; element spacing = lambda/2
Mv      = 32;            % vertical elements
Mh      = 16;            % horizontal elements (N = 512, single polarisation)
BS_H    = 25;            % BS height [m]
UE_H    = 1.5;           % UE height [m]
bore_az = 0;             % global azimuth of array broadside [deg]

beam_name = {'C1','C2','C3','C4','C5','C6','C7','C8','C9','C10','C11','C12','C13','C14'};
ele_tab   = [45.43 36.28 27.8 20.39 14.59 11.14 32 41.55 50.37 12.93 2.29 8.82 15.58 23.35];
azi_tab   = [188.62 184.66 195.4 191.44 222.23 290.9 357.42 358.16 359.61 185.17 24.74 7.02 6.66 5.33];

% ---- angle source: 'sheet' = spreadsheet azimuths, 'photo' = azimuth read from panorama
% Select with BEAM_ANGLE_SOURCE=sheet; default remains photo.
angle_source = getenv('BEAM_ANGLE_SOURCE');
if isempty(angle_source), angle_source = 'photo'; end
assert(any(strcmp(angle_source, {'photo','sheet'})), 'BEAM_ANGLE_SOURCE must be photo or sheet');
px_label = [1497 1392 1225 1137 884 690 268 178 93 1012 787 688 583 483];  % label x-pixel in panorama
F_photo  = 120;      % assumed horizontal FOV of the panorama [deg]
W_photo  = 1602;     % panorama width [px]
xc_photo = 801;      % panorama centre [px] = array broadside (assumed)
if strcmp(angle_source, 'photo')
    azi_tab = (px_label - xc_photo) / W_photo * F_photo;   % right of centre = positive azimuth
end
tag = angle_source;

K = numel(beam_name);
N = Mv*Mh;
c0 = 299792458; lambda = c0/fc; k0 = 2*pi/lambda; d = lambda/2;

%% ---- array geometry: x = 0 plane, y horizontal, z vertical ---------------
[Yg, Zg] = meshgrid(((0:Mh-1) - (Mh-1)/2)*d, ((0:Mv-1) - (Mv-1)/2)*d);
p = [zeros(1,N); Yg(:).'; Zg(:).'];

%% ---- beam weights -------------------------------------------------------
azi_loc = azi_tab - bore_az;
[W, U] = build_beam_codebook(ele_tab, azi_loc, p, fc);

%% ---- sanity: implied horizontal distance if Ele is downtilt -------------
dist_impl = (BS_H - UE_H) ./ tand(ele_tab);

%% ---- array factor and total gain on a 1-degree grid ---------------------
az_g = -180:1:180;   el_g = -90:1:90;
nA = numel(az_g);    nE = numel(el_g);
AF = zeros(K, nE, nA);                 % |w_k^H a(u)|^2, peak = N
Wt = W';                               % K x N
for ie = 1:nE
    el = el_g(ie);
    uu = [cosd(el)*cosd(az_g); cosd(el)*sind(az_g); sind(el)*ones(1,nA)];
    A  = exp(1j*k0*(p.'*uu));          % N x nA
    AF(:, ie, :) = reshape(abs(Wt*A).^2, [K 1 nA]);
end
[AZg, ELg] = meshgrid(az_g, el_g);
AEg = elem_gain_dB(AZg, ELg);          % nE x nA, dBi
TG  = 10*log10(AF) + reshape(AEg, [1 nE nA]);   % total gain, dBi

%% ---- per-beam metrics ---------------------------------------------------
az_req = mod(azi_loc + 180, 360) - 180;
el_req = -ele_tab;
G_peak = 10*log10(N) + elem_gain_dB(az_req, el_req);   % total gain at the commanded direction
bw_az  = zeros(1,K);  bw_el = zeros(1,K);
pk_az  = zeros(1,K);  pk_el = zeros(1,K);  pk_G = zeros(1,K);
cutAz  = zeros(K, nA);
for k = 1:K
    [~, ie0] = min(abs(el_g - el_req(k)));
    [~, ia0] = min(abs(az_g - az_req(k)));
    % -- azimuth 3 dB width at beam elevation (circular), 1 deg grid
    cA = reshape(AF(k, ie0, :), 1, nA);
    cutAz(k,:) = reshape(TG(k, ie0, :), 1, nA);
    thr = N/2;
    [~, i0] = max(cA);  c = circshift(cA, [0, ceil(nA/2) - i0]);  cc = ceil(nA/2);
    l = cc; while l > 1 && c(l-1) >= thr, l = l-1; end
    r = cc; while r < nA && c(r+1) >= thr, r = r+1; end
    bw_az(k) = r - l + 1;
    % -- elevation 3 dB width at beam azimuth (linear)
    cE = reshape(AF(k, :, ia0), 1, nE);
    [~, j0] = max(cE);
    l = j0; while l > 1 && cE(l-1) >= thr, l = l-1; end
    r = j0; while r < nE && cE(r+1) >= thr, r = r+1; end
    bw_el(k) = r - l + 1;
    % -- actual global maximum of the total gain (may be a front-back twin)
    TGk = reshape(TG(k,:,:), nE, nA);
    [tgmax, idx] = max(TGk(:));
    [iep, iap] = ind2sub(size(TGk), idx);
    pk_az(k) = az_g(iap);  pk_el(k) = el_g(iep);  pk_G(k) = tgmax;
end
% mirror twin across the array plane (x -> -x): az -> 180 - az, same array factor
az_mirror = mod(180 - az_req + 180, 360) - 180;

%% ---- inter-beam crosstalk (array factor), dB relative to own peak --------
Gx = abs(Wt * exp(1j*k0*(p.'*U))).^2;  % Gx(k,j) = |w_k^H a_j|^2, diag = N
AE_cmd = elem_gain_dB(az_req, el_req);             % element gain at each commanded direction
XT = 10*log10(Gx / N) + AE_cmd(:).' - AE_cmd(:);   % interference from UE j seen by beam k, relative to own target
XTo = XT;  XTo(1:K+1:end) = -Inf;
[worst_xt, worst_j] = max(XTo, [], 2);
% AF-only crosstalk kept for reference (no element pattern)
XT_af = 10*log10(Gx / N);

%% ---- console report -----------------------------------------------------
% sanity: array factor is identical at az and 180-az (planar array in x=0 plane)
ia = find(az_g == 30); im = find(az_g == 150);
fprintf('\nAF symmetry check  AF(30,el=-20)=%.6f  AF(150,el=-20)=%.6f\n', ...
    AF(1,find(el_g==-20),ia), AF(1,find(el_g==-20),im));
fprintf('\n%-4s %6s %8s %8s %9s %9s %8s %8s %8s %8s %9s\n', 'Beam', 'Ele', 'Azi', 'Dist[m]', 'G@cmd[dBi]', 'Gmax[dBi]', 'AzMax', 'ElMax', 'BW_az', 'BW_el', 'Xtalk[dB]');
for k = 1:K
    fprintf('%-4s %6.2f %8.2f %8.1f %9.2f %9.2f %8.1f %8.1f %8d %8d %9.1f\n', beam_name{k}, ele_tab(k), azi_tab(k), ...
        dist_impl(k), G_peak(k), pk_G(k), pk_az(k), pk_el(k), bw_az(k), bw_el(k), worst_xt(k));
end
fprintf('\nWorst interferer beam per beam (total-gain crosstalk):\n');
for k = 1:K
    fprintf('  %-4s <- %-4s  %6.1f dB   mirror twin of %s: az %.1f\n', beam_name{k}, beam_name{worst_j(k)}, ...
        worst_xt(k), beam_name{k}, az_mirror(k));
end
fprintf('\nArray peak gain (array factor): %.2f dBi, 3GPP element gain: 8 dBi\n', 10*log10(N));
fprintf('Beams with |az|>90 deg (behind the array if broadside is +x): %s\n', ...
    strjoin(beam_name(abs(az_req) > 90), ', '));
fprintf('Beams whose global gain maximum is NOT at the commanded direction (>5 deg): %s\n', ...
    strjoin(beam_name(abs(angle_diff_deg(pk_az, az_req)) > 5 | abs(pk_el - el_req) > 5), ', '));

%% ---- export ----------------------------------------------------------------
save('-v7', fullfile(out, ['beam_model_' tag '.mat']), 'W', 'U', 'p', 'fc', 'N', 'beam_name', 'ele_tab', 'azi_tab', 'XT', 'XT_af');
fid = fopen(fullfile(out, ['beam_summary_' tag '.csv']),'w');
fprintf(fid, 'beam,ele_deg,azi_deg,dist_m,gain_at_cmd_dBi,gmax_dBi,azmax_deg,elmax_deg,bw_az_deg,bw_el_deg,nbr_beam,xtalk_dB,mirror_az_deg\n');
for k = 1:K
    fprintf(fid, '%s,%.2f,%.2f,%.1f,%.2f,%.2f,%.1f,%.1f,%d,%d,%s,%.1f,%.1f\n', beam_name{k}, ele_tab(k), azi_tab(k), ...
        dist_impl(k), G_peak(k), pk_G(k), pk_az(k), pk_el(k), bw_az(k), bw_el(k), beam_name{worst_j(k)}, worst_xt(k), az_mirror(k));
end
fclose(fid);
dlmwrite(fullfile(out, ['beam_cuts_az_' tag '.csv']), cutAz, 'precision', '%.3f');
dlmwrite(fullfile(out, ['beam_crosstalk_' tag '.csv']), XT, 'precision', '%.2f');
printf('\nSaved files with tag %s\n', tag);
