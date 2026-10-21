function report=verify_quadriga1970(qd_path)
%VERIFY_QUADRIGA1970 Test measured power import, physical fields and RX combining.
if nargin>0 && ~isempty(qd_path), addpath(genpath(qd_path)); end
here=fileparts(mfilename('fullpath')); addpath(here);addpath(fullfile(here,'..','common'));
M=load(fullfile(here,'results','ellipse_physical_model.mat'));
files=dir(fullfile(here,'..','..','..','data','measured_patterns','pattern_*.txt'));
assert(numel(files)==78);
[theta,phi]=ndgrid(0:90,0:359);
az=phi(:).'*pi/180; el=(90-theta(:).')*pi/180;
source_error=0;
for i=1:numel(files)
    angles=sscanf(files(i).name,'pattern_Pitch_%d_Azimuth_%d.txt');
    source=dlmread(fullfile(files(i).folder,files(i).name),'\t',6,0);
    for pol={'LHCP','RHCP'}
        a=build_measured_beam1970(angles(1),angles(2),0,pol{1});
        [v,h]=a.interpolate(az,el);
        db=10*log10(abs(v).^2+abs(h).^2);
        expected=source(:,2:end);
        source_error=max(source_error,max(abs(db-expected(:).')));
    end
    if mod(i,13)==0, fprintf('Measured import verified %d/78 (both CP states)\n',i); end
end
assert(source_error<1e-7,'Measured beam power import mismatch');

% Sample directions include strong and weak lobes and every tested scan range.
angles_theta=[0:2:88 1:55]; angles_phi=mod((1:numel(angles_theta))*137,360);
u=[sind(angles_theta).*cosd(angles_phi); ...
   sind(angles_theta).*sind(angles_phi);cosd(angles_theta)];
field_error=0; distance_error=0;
for index=[1 18 31 48 72]
    for pol={'LHCP','RHCP'}
        a=build_physical_array1970(index,[],1,pol{1});
        assert(a.no_elements==M.Nsite && size(a.coupling,1)==M.Nsite);
        assert(size(a.coupling,2)==1,'One link direction has one selected CP beam port');
        [v,h,dist]=a.interpolate(angles_phi*pi/180,(90-angles_theta)*pi/180);
        distance_error=max(distance_error,max(max(abs(dist+M.p_site.'*u))));
        phase=exp(-2j*pi*M.fc/299792458*dist);
        power=abs(a.coupling.'*(v.*phase)).^2+abs(a.coupling.'*(h.*phase)).^2;
        steering=exp(2j*pi*M.fc/299792458*(M.p_site.'*u));
        expected=abs(M.W_receive(:,index)'*steering).^2.*cosd(angles_theta).^M.element_cos_power;
        field_error=max(field_error,norm(power-expected)/max(norm(expected),realmin));
    end
end
assert(field_error<1e-9 && distance_error<1e-10);
% Orthogonality is a property of our ASSUMED CP Jones states, not a measured XPD.
left=build_physical_array1970(1,[],1,'LHCP');
right=build_physical_array1970(1,[],1,'RHCP');
[vl,hl]=left.interpolate(0,pi/4,1);
[vr,hr]=right.interpolate(0,pi/4,1);
cp_inner=abs(conj(vl)*vr+conj(hl)*hr)/(abs(vl)^2+abs(hl)^2);
assert(cp_inner<1e-12);

% Verify actual UPLINK receive combining with a LOS-only channel and exactly
% Nsite RX elements. Reuse one builder/parameter realization; only change Cr.
index=31; a=build_physical_array1970([],[],1,'RHCP');
beam=build_physical_array1970(index,[],1,'RHCP');
s=qd_simulation_parameters; s.center_frequency=M.fc; s.show_progress_bars=0;
s.use_absolute_delays=0;
l=qd_layout(s); l.rx_array=a; l.tx_array=qd_arrayant('omni');
p=M.commands(index,1); az0=M.commands(index,2);
u0=[sind(p)*cosd(az0);sind(p)*sind(az0);cosd(p)];
l.rx_position=[0;0;0]; l.tx_position=1000*u0;
l.set_scenario('LOSonly'); rng(1970,'twister');
b=l.init_builder; b.gen_parameters;
raw=b.get_channels;
assert(size(raw.coeff,1)==M.Nsite && size(raw.coeff,2)==1);
b.rx_array.coupling=beam.coupling;
combined=b.get_channels;
expected=beam.coupling.'*reshape(raw.coeff,M.Nsite,[]);
actual=reshape(combined.coeff,1,[]);
channel_error=norm(actual-expected)/max(norm(expected),realmin);
assert(all(isfinite(actual(:))) && norm(actual)>0 && channel_error<1e-9);
report=struct('measured_files',78,'polarization_states_tested',2, ...
    'measured_power_grid_max_error_db',source_error, ...
    'physical_field_relative_error',field_error,'projection_distance_error_m',distance_error, ...
    'uplink_receive_elements',M.Nsite,'uplink_combined_signal_ports',1, ...
    'uplink_channel_combining_relative_error',channel_error, ...
    'assumed_cp_cross_inner_product_relative',cp_inner, ...
    'absolute_gain_validated',false,'measured_polarization_validated',false, ...
    'measured_complex_phase_validated',false);
save(fullfile(here,'results','quadriga_verification.mat'),'-struct','report','-v7');
fid=fopen(fullfile(here,'results','quadriga_verification.json'),'w');
assert(fid>=0); cleaner=onCleanup(@() fclose(fid));
fprintf(fid,'%s\n',jsonencode(report));
disp(report);
fprintf('PASS: all measured powers, physical array fields, and real uplink RX combining\n');
end
