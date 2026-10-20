function report=verify_allfit_time_domain(qd_path,waveform_file)
%VERIFY_ALLFIT_TIME_DOMAIN Verify directions and the complete 4-slot link.
if nargin<2,waveform_file='';end
here=fileparts(mfilename('fullpath'));addpath(here);addpath(genpath(qd_path));
out=fullfile(here,'results');if ~exist(out,'dir'),mkdir(out);end

% A measured scan node must be preserved exactly before layout orientation.
R=[0 0 1;0 1 0;-1 0 0];
u0=[sind(20)*cosd(30);sind(20)*sind(30);cosd(20)];world=R*u0;
ele=-asind(world(3));azi=atan2d(world(2),world(1));
[~,node,ni]=make_allfit_link_arrays(ele,azi,'RHCP',0);
source=fullfile(here,'..','..','..','data','measured_patterns','pattern_Pitch_20_Azimuth_30.txt');
raw=dlmread(source,'\t',6,0);target=raw(:,2:end);
[theta,phi]=ndgrid(0:90,0:359);
[v,h]=node.interpolate(phi(:).'*pi/180,(90-theta(:).')*pi/180);
nodal_error=max(abs(10*log10(abs(v).^2+abs(h).^2)-target(:).'));
assert(nodal_error<1e-6,'Measured-node power was altered by scan interpolation');
assert(norm(ni.U_world-world)<1e-12);

% Periodic command convention and assumed CP states.
[~,left,li]=make_allfit_link_arrays(ele,azi,'LHCP',0);
[vl,hl]=left.interpolate(0,pi/4);
[vr,hr]=node.interpolate(0,pi/4);
assert(abs((abs(vl)^2+abs(hl)^2)-(abs(vr)^2+abs(hr)^2))<1e-12);
assert(abs(conj(vl)*vr+conj(hl)*hr)<1e-12);
assert(norm(li.Wq-ni.Wq)<1e-12);

% Run unchanged defaults: 4 slots x 14 symbols, 288 subcarriers, FFT=512.
r=iface_allfit_time_domain(struct('qd_path',qd_path,'output_file',waveform_file));
source_cache=fullfile(here,'..','..','..','src','beam_model_photo','beam_model_photo.mat');
S=load(source_cache);
assert(r.fc==1.97e9 && r.N==502 && r.K==14);
assert(isequal(r.ele_tab,S.ele_tab) && isequal(r.azi_tab,S.azi_tab));
assert(r.nSlot==4 && r.nSym==14 && r.nSC==288 && r.nFFT==512 && r.SCS==120e3);
assert(r.nFull==3168 && r.BS_H==25 && r.UE_H==1.5 && r.opts.fd_half==16);
assert(all(r.tot==32256) && all(isfinite(r.ber_meas)) && all(r.ber_meas>=0 & r.ber_meas<=1));
assert(r.consistency_error<1e-10 && r.time_equivalence_error<1e-10);
assert(isequal(find(r.antenna.extrapolated),[1 2 8 9]));
assert(norm(r.antenna.rotation*r.antenna.U_canonical-r.antenna.U_world,'fro')<1e-12);
assert(norm(r.antenna.p_world-r.antenna.rotation*r.antenna.p_canonical,'fro')<1e-12);

% Inspect actual output-beam peaks on the 1-degree canonical pattern grid.
peak_errors=zeros(1,r.K);
for k=1:r.K
    power=r.antenna.power_grid(:,:,k);[~,i]=max(power(:));[ie,ia]=ind2sub(size(power),i);
    el=r.antenna.el_grid_deg(ie);az=r.antenna.az_grid_deg(ia);
    u=[cosd(el)*cosd(az);cosd(el)*sind(az);sind(el)];
    peak_errors(k)=acosd(max(-1,min(1,dot(u,r.antenna.U_canonical(:,k)))));
end
assert(max(peak_errors)<2.5,'Requested directions were clamped or mis-rotated');
relative_correction=zeros(1,r.K);
for u=1:r.K
    base=r.Wq'*r.h_ant{u};
    assert(norm(base+r.h_delta{u}-r.h_beam_cal{u},'fro')/max(norm(r.h_beam_cal{u},'fro'),realmin)<1e-10);
    relative_correction(u)=norm(r.h_delta{u},'fro')/max(norm(base,'fro'),realmin);
end
assert(max(relative_correction)>1e-3,'Calibration path was not exercised');

% Test the retained antenna-noise model, including cross-beam correlations.
rng(1971,'twister');z=(randn(r.N,8192)+1j*randn(r.N,8192))/sqrt(2);
y=r.Wq'*z;covariance=y*y'/size(y,2);expected=r.Wq'*r.Wq;
noise_covariance_error=norm(covariance-expected,'fro')/norm(expected,'fro');
assert(noise_covariance_error<0.08 && min(real(eig(expected)))>-1e-12);

report=struct('fc_hz',r.fc,'receive_elements',r.N,'beams',r.K,'users',r.K, ...
    'slots',r.nSlot,'symbols_per_slot',r.nSym,'subcarriers',r.nSC,'fft',r.nFFT, ...
    'cp_samples',r.nCP,'sample_rate_hz',r.Fs,'bits_per_user',r.tot(1), ...
    'measured_node_max_error_db',nodal_error,'max_peak_direction_error_deg',max(peak_errors), ...
    'time_domain_equivalence_error',r.time_equivalence_error, ...
    'time_frequency_consistency_error',r.consistency_error, ...
    'noise_covariance_relative_error',noise_covariance_error, ...
    'extrapolated_beams',[1 2 8 9],'gain_offset_dbi',r.opts.gain_offset_dbi, ...
    'antenna_phase_is_measured',false,'output_noise_model','unit-norm element combining', ...
    'original_directions_preserved',true,'original_waveform_parameters_preserved',true);
fid=fopen(fullfile(out,'validation.json'),'w');assert(fid>=0);
fprintf(fid,'%s\n',jsonencode(report));fclose(fid);
fid=fopen(fullfile(out,'beam_metrics.csv'),'w');assert(fid>=0);
fprintf(fid,'beam,ele_deg,azi_deg,fit_scan_theta_deg,fit_scan_phi_deg,extrapolated,peak_error_deg,distance_m,mean_sinr_db,ber_pred,ber_meas,error_bits,total_bits\n');
for k=1:r.K
    fprintf(fid,'%s,%.8f,%.8f,%.8f,%.8f,%d,%.8f,%.8f,%.8f,%.8f,%.8f,%d,%d\n', ...
        r.beam_name{k},r.ele_tab(k),r.azi_tab(k),r.antenna.scan_theta_deg(k),r.antenna.scan_phi_deg(k), ...
        r.antenna.extrapolated(k),peak_errors(k),r.dh(k),10*log10(mean(r.sinr(k,:))), ...
        r.ber_pred(k),r.ber_meas(k),r.err(k),r.tot(k));
end
fclose(fid);
save(fullfile(out,'receiver_reference.mat'),'-struct','r','Heff','sinr','ber_pred','ber_meas','Wq','noise_covariance','-v7');
disp(report);fprintf('PASS: directions, calibrated TD transport, default link, and noise model\n');
end
