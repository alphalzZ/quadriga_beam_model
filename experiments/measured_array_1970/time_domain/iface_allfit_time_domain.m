function result = iface_allfit_time_domain(opts)
%IFACE_ALLFIT_TIME_DOMAIN Calibrated 1.97 GHz array in the original OFDM link.
%   result=iface_allfit_time_domain(struct('qd_path','/path/to/QuaDRiGa'));
% All waveform propagation, antenna combining, and residual calibration use
% time samples. The all-fit calibrated model is an output-beam model, not
% measured per-element data: y = Wq' * r_ant + delta_channel * tx.
% The residual FIR makes this equal to the calibrated QuaDRiGa beam response.
if nargin < 1, opts = struct(); end
assert(isstruct(opts) && isscalar(opts), 'opts must be a scalar struct');
defaults = struct('qd_path',getenv('QD_PATH'), 'nSlot',4, 'nSym',14, ...
    'nFFT',512, 'fd_half',16, 'seed',1, 'output_file',fullfile(fileparts(mfilename('fullpath')),'results','iface_allfit_time_domain_out.mat'), ...
    'polarization','RHCP','gain_offset_dbi',0);
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
    candidate = fullfile(here,'..','..','..','QuaDRiGa');
    if exist(candidate,'dir'), addpath(genpath(candidate)); end
end
assert(exist('qd_layout','file') ~= 0, ...
    'QuaDRiGa missing: set QD_PATH, opts.qd_path, or add its sources to the path');
rng(opts.seed,'twister');

%% Sheet fixed-normal Pitch/Azimuth; fitted antenna and carrier.
source_cache=fullfile(here,'..','..','..','src','beam_codebook','results','beam_model_sheet.mat');
assert(exist(source_cache,'file')==2, 'Generate beam_model_sheet.mat with BEAM_ANGLE_SOURCE=sheet first');
S=load(source_cache,'ele_tab','azi_tab','beam_name');
ele=S.ele_tab;azi=S.azi_tab;K=numel(ele);
[physical,calibrated,antenna]=make_allfit_link_arrays(ele,azi,opts.polarization,opts.gain_offset_dbi);
fc=antenna.fc;N=antenna.N;Wq=antenna.Wq;perm=1:N;
BS_H=25;UE_H=1.5;
% Intersect each world-frame beam ray with the UE plane. Pitch is measured
% from the downward array normal, so horizontal range = height gap * tan(Pitch).
ray_distance=(BS_H-UE_H)./(-antenna.U_world(3,:));
ue_positions=bsxfun(@plus,[0;0;BS_H],bsxfun(@times,antenna.U_world,ray_distance));
dh=sqrt(sum(ue_positions(1:2,:).^2,1));
s=qd_simulation_parameters;
s.center_frequency=fc;s.use_absolute_delays=0;s.show_progress_bars=0;
l=qd_layout(s);l.no_tx=1;l.tx_position=[0;0;BS_H];
l.tx_array=physical;l.tx_track.orientation=antenna.orientation;
l.no_rx=K;l.rx_array=qd_arrayant('omni');
for u=1:K
    l.rx_position(:,u)=ue_positions(:,u);
end
l.set_scenario('3GPP_38.901_UMi_NLOS');
% Same scatterer/phase realization for the element and calibrated beam models.
b=l.init_builder;assert(numel(b)==1,'Expected a single common scenario builder');
b.gen_parameters;c_ant=b.get_channels;
% A single builder keeps one TX-array entry per UE/track, even with one BS.
for mobile=1:numel(b.tx_array),b.tx_array(1,mobile)=calibrated;end
c_cal=b.get_channels;
assert(size(c_ant(1).coeff,2)==N && size(c_cal(1).coeff,2)==K);
assert(norm(physical.element_position-antenna.p_canonical,'fro')<1e-12);
% BS -> UE port coefficients are transposed for reciprocal UE -> BS use;
% no conjugation of propagation coefficients is introduced.

%% Sampled impulse responses, including fractional-sample path delays.
nSC = 288; SCS = 120e3; nFull = 3168;
nFFT = opts.nFFT; Fs = nFFT*SCS;
assert(mod(nFFT,2)==0 && nFFT > nSC, 'nFFT must be even and exceed 288');
carrier_index = [-nSC/2:-1 1:nSC/2];        % centered spectrum, DC unused
active = mod(carrier_index,nFFT)+1;         % unshifted FFT bin indices
h_ant=cell(1,K);h_beam_cal=cell(1,K);h_delta=cell(1,K);
delay_origin_s=zeros(1,K);
for u=1:K
    [ca,ta]=paths(c_ant(u),N);
    [cc,tc]=paths(c_cal(u),K);
    delay_origin_s(u)=min([ta(:);tc(:)]);
    h_ant{u}=fractional_fir(ca,max(0,(ta-delay_origin_s(u))*Fs),opts.fd_half);
    h_beam_cal{u}=fractional_fir(cc,max(0,(tc-delay_origin_s(u))*Fs),opts.fd_half);
end
L=max([cellfun(@(h)size(h,2),h_ant),cellfun(@(h)size(h,2),h_beam_cal)]);
nCP=L-1;
assert(nCP<=nFFT,'Delay spread exceeds FFT duration; increase opts.nFFT');
for u=1:K
    h_ant{u}(:,end+1:L)=0;h_beam_cal{u}(:,end+1:L)=0;
    h_delta{u}=h_beam_cal{u}-Wq'*h_ant{u};
end
fprintf('All-data fitted antenna: %.2f GHz, %d RX elements, %d calibrated beams\n',fc/1e9,N,K);
fprintf('Scan extrapolation flags: %s\n',mat2str(antenna.extrapolated));
fprintf('Fs=%.2f MHz, FFT=%d, CP=%d samples\n',Fs/1e6,nFFT,nCP);

%% Perfect receiver CSI from the SAME discretized FIR used by the link.
% This does not apply the channel to the transmitted waveform.
Heff = zeros(K,nSC,K);
for u = 1:K
    H = fft(h_beam_cal{u},nFFT,2);
    Heff(:,:,u) = H(:,active); % already calibrated/combined: no second Wq'
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
consistency_error = 0;
time_equivalence_error = 0;
% Each row selects x[n-tap] for one FIR tap. No circular convolution.
sample_index = bsxfun(@minus,(L:L+nBlock-1), (0:L-1).');
for slot = 1:nSlot
    bits = randi([0 1],K,nSC,nSym,2);
    X = ((1-2*bits(:,:,:,1))+1j*(1-2*bits(:,:,:,2)))/sqrt(2);
    grid = zeros(K,nFFT,nSym);
    grid(:,active,:) = sqrt(Psc)*X;
    tx = ifft(grid,[],2)*sqrt(nFFT);       % unitary OFDM normalization
    tx = reshape(cat(2,tx(:,end-nCP+1:end,:),tx),K,nBlock);
    r_ant = zeros(N,nBlock);              % 502 physical antenna time streams
    y_correction = zeros(K,nBlock);
    y_direct = zeros(K,nBlock);
    for u = 1:K
        extended = [history(u,:) tx(u,:)];
        delayed_samples = extended(sample_index); % L taps x nBlock samples
        % r_a[n] += sum_t h_a,u[t] * x_u[n-t], for all antennas at once.
        r_ant = r_ant + h_ant{u}*delayed_samples;
        y_correction = y_correction + h_delta{u}*delayed_samples;
        y_direct = y_direct + h_beam_cal{u}*delayed_samples;
        history(u,:) = extended(end-L+2:end);
    end
    y_clean=Wq'*r_ant+y_correction;
    td_error=norm(y_clean(:)-y_direct(:))/max(norm(y_direct(:)),realmin);
    time_equivalence_error=max(time_equivalence_error,td_error);
    assert(td_error<1e-10,'Element propagation plus time residual must equal calibrated FIR');
    if slot == 1 || slot == nSlot
        % Check CP/stream continuity against the calibrated receiver CSI.
        Yclean = demodulate(y_clean,nFFT,nCP,nSym,active);
        Yref = zeros(K,nSC,nSym);
        for u = 1:K
            Yref = Yref + bsxfun(@times,g(:,:,u),reshape(X(u,:,:),1,nSC,nSym));
        end
        fd_error = norm(Yclean(:)-Yref(:))/max(norm(Yref(:)),realmin);
        consistency_error=max(consistency_error,fd_error);
        assert(consistency_error < 1e-10, 'Time-domain channel consistency check failed');
        fprintf('Noiseless time/frequency reference relative error: %.3g\n',consistency_error);
    end
    % Unit-norm beams preserve variance N0. Antenna-domain white noise also
    % preserves cross-beam covariance N0*(Wq''*Wq), unlike independent beams.
    % With unitary IFFT/FFT, normalized sample noise variance is also N0.
    r_ant = r_ant + sqrt(N0/2)*(randn(N,nBlock)+1j*randn(N,nBlock));
    y_beam = Wq'*r_ant + y_correction;     % element combining + calibrated TIME FIR
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
result = struct('h_ant',{h_ant},'h_beam_cal',{h_beam_cal},'h_delta',{h_delta}, ...
    'source_cache',source_cache,'antenna',antenna,'ele_tab',ele,'azi_tab',azi,'beam_name',{S.beam_name}, ...
    'BS_H',BS_H,'UE_H',UE_H,'ue_positions',ue_positions,'dh',dh,'N',N,'K',K,'nFull',nFull, ...
    'time_equivalence_error',time_equivalence_error,'noise_covariance',N0*(Wq'*Wq), ...
    'Wq',Wq,'perm',perm,'Heff',Heff, ...
    'g',g,'fc',fc,'Fs',Fs,'SCS',SCS,'nFFT',nFFT,'nCP',nCP,'nSC',nSC, ...
    'nSym',nSym,'nSlot',nSlot,'active',active,'carrier_index',carrier_index, ...
    'Psc',Psc,'N0',N0,'sinr',sinr,'ber_pred',ber_pred,'ber_meas',ber_meas, ...
    'err',err,'tot',tot,'delay_origin_s',delay_origin_s, ...
    'fractional_delay_latency_samples',opts.fd_half,'opts',opts, ...
    'consistency_error',consistency_error,'tx_last',tx,'y_beam_last',y_beam);
if ~isempty(opts.output_file)
    output_parent=fileparts(opts.output_file);
    if ~isempty(output_parent) && ~exist(output_parent,'dir'),mkdir(output_parent);end
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

function [coeff,tau]=paths(channel,N)
coeff=reshape(channel.coeff(1,:,:,1),N,[]);
if channel.individual_delays
    tau=reshape(channel.delay(1,:,:,1),N,[]);
else
    tau=repmat(channel.delay(:,1).',N,1);
end
assert(isequal(size(coeff),size(tau)) && all(isfinite(coeff(:))) && all(isfinite(tau(:))));
end
