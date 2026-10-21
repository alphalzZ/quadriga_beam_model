% IFACE_DEMO  QuaDRiGa channel -> beam-domain effective channel -> multi-user
%             uplink symbol transport -> per-beam MMSE -> BER, compared with SINR prediction
here=fileparts(mfilename('fullpath'));
QD_PATH=getenv('QD_PATH');
if isempty(QD_PATH),QD_PATH=fullfile(here,'..','..','QuaDRiGa');end
out=fullfile(here,'results');if ~exist(out,'dir'),mkdir(out);end
addpath(genpath(QD_PATH));
more off;
rand('state', 1); randn('state', 1);

%% ---- beams (photo-corrected) -------------------------------------------
S   = load(fullfile(here,'..','beam_codebook','results','beam_model_photo.mat'));             % beam weights from run_beam_model (angle_source='photo')
fc  = S.fc;  W = S.W;  pm = S.p;  N = S.N;          % W: 512 x 14, unit-norm columns
ele = S.ele_tab;  azi = S.azi_tab;  K = numel(ele);
BS_H = 25;  UE_H = 1.5;
dh = (BS_H - UE_H) ./ tand(ele);                    % horizontal distance implied by downtilt

%% ---- OFDM numerology ----------------------------------------------------
SCS  = 120e3;                    % subcarrier spacing
nSC  = 288;                      % 24 PRB simulated (demo; full band = 3168)
nSym = 14;                       % symbols per slot
nFull = 3168;                    % PSD reference: 264 PRB full band
c0 = 299792458;

%% ---- QuaDRiGa layout: UEs placed on their beam directions --------------
s = qd_simulation_parameters;
s.center_frequency    = fc;
s.use_absolute_delays = 0;       % relative delays (UE timing advance removed)
s.show_progress_bars  = 0;
l = qd_layout(s);
l.no_tx = 1;  l.tx_position = [0; 0; BS_H];
l.tx_array = qd_arrayant('3gpp-3d', 32, 16, fc, 1);
l.no_rx = K;  l.rx_array = qd_arrayant('omni');
for k = 1:K
    l.rx_position(:, k) = [dh(k)*cosd(azi(k)); dh(k)*sind(azi(k)); UE_H];
end
l.set_scenario('3GPP_38.901_UMi_NLOS');
c = l.get_channels;

%% ---- port ordering check: QuaDRiGa element order vs. beam weights --------
pq = l.tx_array.element_position;
perm = zeros(1, 512);
for i = 1:512
    d2 = sum((pm - repmat(pq(:,i), 1, 512)).^2, 1);
    [dmin, perm(i)] = min(d2);
    if dmin > 1e-12, error('element geometry mismatch at port %d', i); end
end
fprintf('port permutation check: unique = %d (1 = identity ordering)\n', numel(unique(perm)));
Wq = W(perm, :);                 % weights aligned with QuaDRiGa port order

%% ---- frequency response on the 120 kHz grid -----------------------------
BWd = nSC * SCS;                 % fr() carrier spacing = BW/nSC = SCS
Hf = zeros(512, nSC, K);
for k = 1:K
    Hf(:,:,k) = reshape(fr(c(k), BWd, nSC), 512, nSC);
end
fprintf('carrier spacing check: %.1f kHz\n', BWd/nSC/1e3);

%% ---- beam-domain effective channel -------------------------------------
Heff = zeros(K, nSC, K);         % (beam, subcarrier, UE)
for u = 1:K
    Heff(:,:,u) = Wq' * Hf(:,:,u);
end

%% ---- powers -------------------------------------------------------------
kB = 1.38e-23;  T = 290;  NF = 10^(7/10);
N0 = kB * T * NF * SCS;                         % noise per RE, per beam (unit-norm W)
Psc = 10^((23-30)/10) / nFull;                  % UE PSD of a full-band allocation
g = sqrt(Psc) * Heff;                           % scaled effective channel

%% ---- analytic prediction per UE (Gaussian interference) ----------------
sinr = zeros(K, nSC);  sinr_snr = zeros(K, nSC);
for u = 1:K
    gu = reshape(g(u,:,u), nSC, 1);
    Iu = reshape(sum(abs(g(u,:,setdiff(1:K,u))).^2, 3), nSC, 1) + N0;
    sinr(u,:)     = (abs(gu).^2 ./ Iu).';
    sinr_snr(u,:) = (abs(gu).^2 / N0).';
end
Qf = @(x) 0.5 * erfc(x / sqrt(2));
ber_pred     = mean(Qf(sqrt(sinr)), 2);
ber_pred_snr = mean(Qf(sqrt(sinr_snr)), 2);

%% ---- waveform-level transmission and receiver ---------------------------
nSlot = 40;
err = zeros(K,1);  tot = zeros(K,1);
for slot = 1:nSlot
    bits = randi([0 1], K, nSC, nSym, 2);
    X = zeros(K, nSC, nSym);
    for u = 1:K
        X(u,:,:) = ((1 - 2*reshape(bits(u,:,:,1), 1, nSC, nSym)) + ...
                    1j*(1 - 2*reshape(bits(u,:,:,2), 1, nSC, nSym))) / sqrt(2);
    end
    Y = zeros(K, nSC, nSym);
    for u = 1:K
        Y = Y + reshape(g(:,:,u), K, nSC, 1) .* reshape(X(u,:,:), 1, nSC, nSym);
    end
    Y = Y + sqrt(N0/2) * (randn(K, nSC, nSym) + 1j*randn(K, nSC, nSym));
    for u = 1:K
        y  = reshape(Y(u,:,:), nSC, nSym);          % beam u
        gu = reshape(g(u,:,u), nSC, 1);             % nSC x 1
        Iu = reshape(sum(abs(g(u,:,setdiff(1:K,u))).^2, 3), nSC, 1) + N0;
        xh = (conj(gu) ./ (abs(gu).^2 + Iu)) .* y; % per-RE MMSE with interference as noise
        bh_I = real(xh) < 0;  bh_Q = imag(xh) < 0; % Gray QPSK hard decision
        b1 = reshape(bits(u,:,:,1), nSC, nSym);
        b2 = reshape(bits(u,:,:,2), nSC, nSym);
        err(u) = err(u) + sum(sum(bh_I ~= b1)) + sum(sum(bh_Q ~= b2));
        tot(u) = tot(u) + 2*nSC*nSym;
    end
end
ber_meas = err ./ tot;

%% ---- report -------------------------------------------------------------
fprintf('\nUE  dist[m]  meanSINR[dB]  BER_pred(I+N)  BER_pred(N only)  BER_meas\n');
for u = 1:K
    fprintf('%2d  %7.1f  %10.1f  %13.4f  %16.4f  %9.4f\n', u, dh(u), ...
        10*log10(mean(sinr(u,:))), ber_pred(u), ber_pred_snr(u), ber_meas(u));
end
save('-v7', fullfile(out,'iface_out.mat'), 'Heff', 'g', 'Psc', 'N0', 'nSC', 'nSym', 'SCS', 'ber_pred', 'ber_meas', 'sinr');
printf('saved iface_out.mat\n');
