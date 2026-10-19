function result = iface_demo_time_domain(opts)
% IFACE_DEMO_TIME_DOMAIN  Sample-domain multipath and receive beamforming.
%   result = iface_demo_time_domain();
%   result = iface_demo_time_domain(struct('qd_path','/path/to/QuaDRiGa', ...
%       'nSlot',1,'nSym',2,'output_file','/tmp/iface_td.mat'));
%
% Signal path: QPSK -> IFFT -> CP -> per-antenna FIR convolution -> antenna
% noise -> Wq' * received samples -> remove CP -> FFT -> scalar MMSE -> BER.
% No frequency-domain channel multiplication is used to transport samples.
% The static, frequency-flat beam weights act on every time sample. The FFT
% of the FIR is used ONLY for receiver CSI and a noiseless consistency check.
% Requires GNU Octave or MATLAB and QuaDRiGa 2.8.1 (no signal toolbox needed).
if nargin < 1, opts = struct(); end
assert(isstruct(opts) && isscalar(opts), 'opts must be a scalar struct');
defaults = struct('qd_path',getenv('QD_PATH'), 'nSlot',4, 'nSym',14, ...
    'nFFT',512, 'fd_half',16, 'seed',1, 'output_file','iface_time_domain_out.mat');
names = fieldnames(opts);
for i = 1:numel(names)
    assert(isfield(defaults,names{i}), ['Unknown option: ' names{i}]);
    defaults.(names{i}) = opts.(names{i});
end
opts = defaults;
for name = {'nSlot','nSym','nFFT','fd_half','seed'}
    v = opts.(name{1});
    assert(isnumeric(v) && isscalar(v) && isfinite(v) && v == fix(v) && v >= 1, ...
        [name{1} ' must be a positive integer']);
end
here = fileparts(mfilename('fullpath'));
if ~isempty(opts.qd_path), addpath(genpath(opts.qd_path)); end
if exist('qd_layout','file') == 0
    candidate = fullfile(here,'..','..','QuaDRiGa');
    if exist(candidate,'dir'), addpath(genpath(candidate)); end
end
assert(exist('qd_layout','file') ~= 0, ...
    'QuaDRiGa missing: set QD_PATH, opts.qd_path, or add its sources to the path');
rng(opts.seed,'twister');

%% Array and reciprocal, static channel (one snapshot, no Doppler).
S = load(fullfile(here,'beam_model_photo.mat'));
fc = 3.5e9; Mv = 32; Mh = 16; N = Mv*Mh; K = numel(S.ele_tab);
assert(isscalar(S.fc) && abs(S.fc-fc) < 1 && S.N == N, ...
    'Beam cache must use 3.5 GHz and 512 ports; rerun run_beam_model.m');
assert(isequal(size(S.W),[N K]) && isequal(size(S.p),[3 N]));
BS_H = 25; UE_H = 1.5;
dh = (BS_H-UE_H) ./ tand(S.ele_tab);
assert(numel(unique(S.p(3,:))) == Mv && numel(unique(S.p(2,:))) == Mh, ...
    'Beam cache must use vertical 32 / horizontal 16; rerun run_beam_model.m');
s = qd_simulation_parameters;
s.center_frequency = fc; s.use_absolute_delays = 0; s.show_progress_bars = 0;
l = qd_layout(s);
l.no_tx = 1; l.tx_position = [0;0;BS_H];
% QuaDRiGa takes vertical count FIRST, matching run_beam_model.m.
l.tx_array = qd_arrayant('3gpp-3d',Mv,Mh,fc,1);
l.no_rx = K; l.rx_array = qd_arrayant('omni');
for u = 1:K
    l.rx_position(:,u) = [dh(u)*cosd(S.azi_tab(u)); ...
        dh(u)*sind(S.azi_tab(u)); UE_H];
end
l.set_scenario('3GPP_38.901_UMi_NLOS');
pq = l.tx_array.element_position;
assert(isequal(size(pq),size(S.p)), 'Array port count mismatch');
perm = zeros(1,N);
for a = 1:N
    [dmin,perm(a)] = min(sum((S.p-pq(:,a)).^2,1));
    assert(dmin < 1e-12, 'Array element geometry mismatch');
end
assert(numel(unique(perm)) == N, 'Port mapping must be bijective');
Wq = S.W(perm,:);                          % N antennas x K beams
assert(max(abs(sum(abs(Wq).^2,1)-1)) < 1e-10, 'Weights must have unit norm');
c = l.get_channels;
% QuaDRiGa generates BS -> UE coefficients. Reciprocity gives UE -> BS
% coefficients by TRANSPOSE of the port indices, without conjugation.

%% Sampled impulse responses, including fractional-sample path delays.
nSC = 288; SCS = 120e3; nFull = 3168;
nFFT = opts.nFFT; Fs = nFFT*SCS;
assert(mod(nFFT,2)==0 && nFFT > nSC, 'nFFT must be even and exceed 288');
carrier_index = [-nSC/2:-1 1:nSC/2];        % centered spectrum, DC unused
active = mod(carrier_index,nFFT)+1;         % unshifted FFT bin indices
h_ant = cell(1,K); delay_origin_s = zeros(1,K);
for u = 1:K
    coeff = reshape(c(u).coeff(1,:,:,1),N,[]);
    if c(u).individual_delays
        tau = reshape(c(u).delay(1,:,:,1),N,[]);
    else
        tau = repmat(c(u).delay(:,1).',N,1);
    end
    assert(isequal(size(coeff),size(tau)), 'Unexpected QuaDRiGa delay dimensions');
    assert(all(isfinite(coeff(:))) && all(isfinite(tau(:))));
    % Ideal per-UE timing advance: remove one common delay, preserving path
    % excess delays AND relative antenna delays. fd_half adds causal latency.
    delay_origin_s(u) = min(tau(:));
    delays = max(0,(tau-delay_origin_s(u))*Fs);
    h_ant{u} = fractional_fir(coeff,delays,opts.fd_half);
end
L = max(cellfun(@(h) size(h,2),h_ant));
nCP = L-1;                               % covers FIR tail incl. interpolation
assert(nCP <= nFFT, 'Delay spread exceeds FFT duration; increase opts.nFFT');
for u = 1:K, h_ant{u}(:,end+1:L) = 0; end
fprintf('Array vertical=%d horizontal=%d, fc=%.2f GHz, %d beams\n',Mv,Mh,fc/1e9,K);
fprintf('Fs=%.2f MHz, FFT=%d, CP=%d samples\n',Fs/1e6,nFFT,nCP);

%% Perfect receiver CSI from the SAME discretized FIR used by the link.
% This does not apply the channel to the transmitted waveform.
Heff = zeros(K,nSC,K);
for u = 1:K
    H = fft(h_ant{u},nFFT,2);
    Heff(:,:,u) = Wq' * H(:,active);
end
kB = 1.38e-23; T = 290; NF = 10^(7/10);
N0 = kB*T*NF*SCS;                         % complex noise variance per RE
Psc = 10^((23-30)/10)/nFull;
g = sqrt(Psc)*Heff;
sinr = zeros(K,nSC); equalizer = zeros(K,nSC);
for u = 1:K
    gu = g(u,:,u);
    interference = sum(abs(g(u,:,setdiff(1:K,u))).^2,3)+N0;
    sinr(u,:) = abs(gu).^2 ./ interference;
    equalizer(u,:) = conj(gu) ./ (abs(gu).^2+interference);
end
ber_pred = mean(0.5*erfc(sqrt(sinr)/sqrt(2)),2);

%% Actual waveform transport: all channel filtering and combining in time.
nSym = opts.nSym; nSlot = opts.nSlot;
nBlock = (nFFT+nCP)*nSym;
history = zeros(K,L-1);                   % continuity across slot boundaries
err = zeros(K,1); tot = zeros(K,1);
consistency_error = NaN;
% Each row selects x[n-tap] for one FIR tap. No circular convolution.
sample_index = bsxfun(@minus,(L:L+nBlock-1), (0:L-1).');
for slot = 1:nSlot
    bits = randi([0 1],K,nSC,nSym,2);
    X = ((1-2*bits(:,:,:,1))+1j*(1-2*bits(:,:,:,2)))/sqrt(2);
    grid = zeros(K,nFFT,nSym);
    grid(:,active,:) = sqrt(Psc)*X;
    tx = ifft(grid,[],2)*sqrt(nFFT);       % unitary OFDM normalization
    tx = reshape(cat(2,tx(:,end-nCP+1:end,:),tx),K,nBlock);
    r_ant = zeros(N,nBlock);              % received time samples BEFORE beams
    for u = 1:K
        extended = [history(u,:) tx(u,:)];
        delayed_samples = extended(sample_index); % L taps x nBlock samples
        % r_a[n] += sum_t h_a,u[t] * x_u[n-t], for all antennas at once.
        r_ant = r_ant + h_ant{u}*delayed_samples;
        history(u,:) = extended(end-L+2:end);
    end
    if slot == 1
        % Verify CP + linear filtering + sample-domain combining against CSI.
        Yclean = demodulate(Wq'*r_ant,nFFT,nCP,nSym,active);
        Yref = zeros(K,nSC,nSym);
        for u = 1:K
            Yref = Yref + bsxfun(@times,g(:,:,u),reshape(X(u,:,:),1,nSC,nSym));
        end
        consistency_error = norm(Yclean(:)-Yref(:))/max(norm(Yref(:)),realmin);
        assert(consistency_error < 1e-10, 'Time-domain channel consistency check failed');
        fprintf('Noiseless time/frequency reference relative error: %.3g\n',consistency_error);
    end
    % Unit-norm beams preserve variance N0. Antenna-domain white noise also
    % preserves cross-beam covariance N0*(Wq''*Wq), unlike independent beams.
    % With unitary IFFT/FFT, normalized sample noise variance is also N0.
    r_ant = r_ant + sqrt(N0/2)*(randn(N,nBlock)+1j*randn(N,nBlock));
    y_beam = Wq'*r_ant;                   % K beams x nBlock TIME samples
    Y = demodulate(y_beam,nFFT,nCP,nSym,active);
    for u = 1:K
        xhat = bsxfun(@times,equalizer(u,:).',reshape(Y(u,:,:),nSC,nSym));
        bI = reshape(bits(u,:,:,1),nSC,nSym);
        bQ = reshape(bits(u,:,:,2),nSC,nSym);
        err(u) = err(u)+sum(sum((real(xhat)<0)~=bI))+sum(sum((imag(xhat)<0)~=bQ));
        tot(u) = tot(u)+2*nSC*nSym;
    end
    fprintf('Slot %d/%d complete\n',slot,nSlot);
end
ber_meas = err./tot;
fprintf('\nUE  dist[m]  meanSINR[dB]  BER_pred(I+N)  BER_meas\n');
for u = 1:K
    fprintf('%2d %8.1f %13.2f %14.5f %10.5f\n',u,dh(u), ...
        10*log10(mean(sinr(u,:))),ber_pred(u),ber_meas(u));
end
% Keep the last slot waveform for inspection; all slots share this static FIR.
result = struct('h_ant',{h_ant},'Wq',Wq,'perm',perm,'Heff',Heff, ...
    'g',g,'fc',fc,'Fs',Fs,'SCS',SCS,'nFFT',nFFT,'nCP',nCP,'nSC',nSC, ...
    'nSym',nSym,'nSlot',nSlot,'active',active,'carrier_index',carrier_index, ...
    'Psc',Psc,'N0',N0,'sinr',sinr,'ber_pred',ber_pred,'ber_meas',ber_meas, ...
    'err',err,'tot',tot,'delay_origin_s',delay_origin_s, ...
    'fractional_delay_latency_samples',opts.fd_half,'opts',opts, ...
    'consistency_error',consistency_error,'tx_last',tx,'y_beam_last',y_beam);
if ~isempty(opts.output_file)
    save(opts.output_file,'-struct','result','-v7');
    fprintf('Saved %s\n',opts.output_file);
end
end

function h = fractional_fir(coeff,delays,half)
% Windowed-sinc fractional delay, 2*half+1 taps per path. Explicit causal
% latency of half samples avoids clipping the earliest path's sinc kernel.
% Normalize each kernel to unit DC gain; this is a finite-band approximation
% to continuous delays, not nearest-sample rounding or an exact delay model.
[N,P] = size(coeff);
L = max(floor(delays(:)))+2*half+1;
h = zeros(N,L);
t = -half:half;
win = 0.5+0.5*cos(pi*t/half);
for a = 1:N
    for p = 1:P
        d = delays(a,p); integer = floor(d); frac = d-integer;
        z = t-frac;
        kernel = ones(size(z)); nz = abs(z)>1e-14;
        kernel(nz) = sin(pi*z(nz))./(pi*z(nz));
        kernel = kernel.*win; kernel = kernel/sum(kernel);
        ix = integer+(1:2*half+1);
        h(a,ix) = h(a,ix)+coeff(a,p)*kernel;
    end
end
end

function Y = demodulate(samples,nFFT,nCP,nSym,active)
K = size(samples,1);
symbols = reshape(samples,K,nFFT+nCP,nSym);
grid = fft(symbols(:,nCP+1:end,:),[],2)/sqrt(nFFT);
Y = grid(:,active,:);
end
