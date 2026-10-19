function report = verify_beam_model_triangular_cp(qd_path)
%VERIFY_BEAM_MODEL_TRIANGULAR_CP Geometry, CP, cache and QuaDRiGa integration.
%   report = verify_beam_model_triangular_cp('/path/to/QuaDRiGa');
if nargin>0 && ~isempty(qd_path), addpath(genpath(qd_path)); end
assert(exist('qd_arrayant','file')~=0,'QuaDRiGa is required for this verification');
base = struct('output_dir','','grid_step_deg',5);
lh = base; lh.polarization='LHCP';
rh = base; rh.polarization='RHCP';
du = base; du.polarization='DUAL';
L = run_beam_model_triangular_cp(lh);
R = run_beam_model_triangular_cp(rh);
D = run_beam_model_triangular_cp(du);
assert(isequal(size(R.W),[512 14]) && isequal(size(D.W),[1024 28]));
assert(norm(mean(R.p_site,2))<1e-14);
d = 299792458/R.fc/2;
p = reshape(R.p_site,3,32,16);
assert(max(abs(diff(squeeze(p(3,:,1)))-sqrt(3)/2*d))<1e-12);
assert(max(abs(abs(diff(squeeze(p(2,:,1))))-d/2))<1e-12);
assert(max(abs(diff(squeeze(p(2,1,:)))-d))<1e-12);
assert(abs(R.nearest_spacing_m-d)<1e-12);
assert(norm(L.W-R.W,'fro')<1e-12, 'CP sign must not alter spatial weights');
assert(norm(L.XT-R.XT,'fro')<1e-10, 'Matched LHCP/RHCP patterns must agree');
% No undesired +3 dB from counting each CP field component as full power.
expected = 10*log10(512)+elem_gain_dB(R.azi_tab,-R.ele_tab);
assert(max(abs(R.gain_cmd-expected))<1e-12);
J = D.port_jones(:,[1 513]);
assert(norm(J'*J-eye(2),'fro')<1e-12);
assert(max(max(D.Gx(1:14,15:28)))<1e-20);
assert(max(max(D.Gx(15:28,1:14)))<1e-20);
assert(max(abs(diag(D.Gx)-512))<1e-8);

% Rectangular compatibility with the independently generated legacy cache.
rect = rh; rect.lattice='rectangular';
Q = run_beam_model_triangular_cp(rect);
old = load(fullfile(fileparts(mfilename('fullpath')),'beam_model_photo.mat'));
assert(old.fc==Q.fc && norm(old.p-Q.p,'fro')<1e-12);
assert(norm(old.W-Q.W,'fro')<1e-10);
bad = base; bad.lattice='not-a-lattice';
caught = false;
try, run_beam_model_triangular_cp(bad); catch, caught=true; end
assert(caught, 'Unsupported lattice must fail');

% Construct the full-size QuaDRiGa adapter and sample non-broadside fields.
a = make_cp_arrayant(D,5);
az = [-45 0 35]*pi/180; el = [-20 0 15]*pi/180;
[V,H,dist] = a.interpolate(az,el);
amp = sqrt(10.^(elem_gain_dB(az*180/pi,el*180/pi)/10));
assert(isequal(size(V),[1024 3]));
assert(norm(V-D.port_jones(1,:).'*amp,'fro')/norm(V,'fro')<1e-10);
assert(norm(H-D.port_jones(2,:).'*amp,'fro')/norm(H,'fro')<1e-10);
U = [cos(el).*cos(az);cos(el).*sin(az);sin(el)];
assert(norm(dist+D.p.'*U,'fro')<1e-10);
assert(max(max(abs(abs(V).^2+abs(H).^2-ones(1024,1)*amp.^2)))<1e-10);
% Cross-check CP handedness at broadside against official H/V element
% generation and its documented coupling convention, rather than labels only.
ref = qd_arrayant('3gpp-3d',1,1,R.fc,2);
ref.coupling = [1 1;1j -1j]/sqrt(2);
ref.combine_pattern;
[vr,hr] = ref.interpolate(0,0);
[va,ha] = a.interpolate(0,0,[1 513]);
assert(norm([va ha]-[vr hr],'fro')<1e-10);

% A real QuaDRiGa channel with a nonrectangular, dual-CP antenna. Use a small
% aperture to keep the smoke test fast; the full 1024-port fields were checked above.
small = du; small.Mv=4; small.Mh=3;
M = run_beam_model_triangular_cp(small);
s = qd_simulation_parameters; s.center_frequency=M.fc; s.show_progress_bars=0;
l = qd_layout(s); l.tx_array=make_cp_arrayant(M,5);
l.rx_array=qd_arrayant('lhcp-rhcp-dipole');
l.tx_position=[0;0;25]; l.rx_position=[100;20;1.5];
l.set_scenario('3GPP_38.901_UMi_LOS');
rng(7,'twister'); c=l.get_channels;
assert(size(c.coeff,1)==2 && size(c.coeff,2)==24);
assert(all(isfinite(c.coeff(:))) && norm(c.coeff(:))>0);

% Export/load round trip in an isolated directory, removed after verification.
out = tempname; mkdir(out); cleanup = onCleanup(@() rmdir(out,'s'));
small.output_dir=out; M=run_beam_model_triangular_cp(small);
cache = load(fullfile(out,'beam_model_triangular_dual.mat'));
assert(isequal(cache.W,M.W) && isequal(cache.port_jones,M.port_jones));
assert(exist(fullfile(out,'beam_model_triangular_dual_summary.csv'),'file')==2);
report = struct('geometry',true,'polarization',true,'rectangular_regression',true, ...
    'quadriga_fields',true,'quadriga_channel',true,'cache_round_trip',true, ...
    'single_ports',512,'dual_ports',1024,'single_beams',14,'dual_beams',28);
disp(report);
fprintf('PASS: geometry, CP isolation/normalization, rectangular regression,\n');
fprintf('      QuaDRiGa fields/channel generation, and cache round trip.\n');
end
