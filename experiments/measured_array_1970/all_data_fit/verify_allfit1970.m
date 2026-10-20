function report=verify_allfit1970(qd_path)
%VERIFY_ALLFIT1970 Check all calibrated beam outputs against all source samples.
if nargin>0 && ~isempty(qd_path),addpath(genpath(qd_path));end
here=fileparts(mfilename('fullpath'));addpath(here);
root=fullfile(here,'..','..','..','data','measured_patterns');
files=dir(fullfile(root,'pattern_*.txt'));assert(numel(files)==78);
[theta,phi]=ndgrid(0:90,0:359);
az=phi(:).'*pi/180;el=(90-theta(:).')*pi/180;
max_error=0;sum_square=0;count=0;
per_file=struct([]);
for i=1:numel(files)
    angles=sscanf(files(i).name,'pattern_Pitch_%d_Azimuth_%d.txt');
    raw=dlmread(fullfile(root,files(i).name),'\t',6,0);target=raw(:,2:end);
    file_max=0;
    for pol={'LHCP','RHCP'}
        a=build_allfit_beam1970(angles(1),angles(2),pol{1});
        [v,h]=a.interpolate(az,el);
        gain=10*log10(abs(v).^2+abs(h).^2);
        delta=gain-target(:).';
        file_max=max(file_max,max(abs(delta)));
        max_error=max(max_error,max(abs(delta)));sum_square=sum_square+sum(delta.^2);count=count+numel(delta);
    end
    per_file(i).filename=files(i).name;per_file(i).max_error_db=file_max;
    if mod(i,13)==0,fprintf('Calibrated QuaDRiGa outputs checked: %d/78\n',i);end
end
assert(max_error<1e-7,'Calibrated output does not reproduce source powers');
% Verify a real channel can consume the calibrated pattern. This is an
% interface/finite-coefficient check, not a measured phase/channel comparison.
s=qd_simulation_parameters;s.center_frequency=1.97e9;s.show_progress_bars=0;
l=qd_layout(s);l.rx_array=build_allfit_beam1970(20,30);l.tx_array=qd_arrayant('omni');
l.rx_position=[0;0;0];
l.tx_position=1000*[sind(20)*cosd(30);sind(20)*sind(30);cosd(20)];
l.set_scenario('LOSonly');rng(1970,'twister');c=l.get_channels;
assert(size(c.coeff,1)==1 && size(c.coeff,2)==1 && norm(c.coeff(:))>0 && all(isfinite(c.coeff(:))));
report=struct('fit_files',78,'holdout_files',0,'source_samples',78*91*360, ...
    'cp_states_tested',2,'total_power_comparisons',count,'quadriga_rmse_db',sqrt(sum_square/count), ...
    'quadriga_max_error_db',max_error,'channel_interface_passed',true, ...
    'measured_phase_validated',false,'measured_polarization_validated',false);
fid=fopen(fullfile(here,'results','quadriga_verification.json'),'w');assert(fid>=0);
cleanup=onCleanup(@() fclose(fid));fprintf(fid,'%s\n',jsonencode(report));
save(fullfile(here,'results','quadriga_verification.mat'),'-struct','report','-v7');
disp(report);fprintf('PASS: calibrated simulation powers match all measured grids\n');
end
