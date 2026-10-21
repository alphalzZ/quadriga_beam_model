function [physical,calibrated,info]=make_allfit_link_arrays(ele_deg,azi_deg,polarization,gain_offset_dbi)
%MAKE_ALLFIT_LINK_ARRAYS 502 element baseline and 14 calibrated beam outputs.
% Canonical fitted broadside is +z; assume a nadir-facing array in the world.
% Bank pi maps +z to -z: R=diag([1 -1 -1]), orientation [pi;0;0].
% Spreadsheet ele is measured off-normal Pitch: theta=ele, phi=azi.
% Construct in this fixed canonical frame, then rotate into world coordinates.
% Measured scan commands are interpolated in LINEAR POWER after transporting
% each corner pattern by the minimal 3D rotation from its command to target.
% Outside the measured 0..50 degree scan range, transport the boundary pattern;
% do not clamp the requested direction or present extrapolation as measured.
if nargin<3,polarization='RHCP';end
if nargin<4,gain_offset_dbi=0;end
assert(any(strcmpi(polarization,{'LHCP','RHCP'})));
ele_deg=reshape(ele_deg,1,[]);azi_deg=reshape(azi_deg,1,[]);
assert(numel(ele_deg)==numel(azi_deg) && all(isfinite([ele_deg azi_deg])));
assert(isscalar(gain_offset_dbi) && isfinite(gain_offset_dbi));
here=fileparts(mfilename('fullpath'));parent=fileparts(here);
addpath(fullfile(parent,'common'));addpath(fullfile(parent,'all_data_fit'));
model_path=fullfile(parent,'all_data_fit','results','physical_baseline.mat');
M=load(model_path);N=M.Nsite;
assert(M.fc==1.97e9 && N==502);
physical=build_physical_array1970([],model_path,1,polarization);
R=diag([1 -1 -1]);
assert(all(ele_deg>=0 & ele_deg<90), 'Pitch must lie in the front hemisphere [0,90)');
theta=ele_deg;phi=mod(azi_deg,360);
U=[sind(theta).*cosd(phi);sind(theta).*sind(phi);cosd(theta)];
U_world=R*U;
K=numel(theta);
W=bsxfun(@times,M.taper(:),exp(2j*pi*M.fc/299792458*(M.p_site.'*U)))/norm(M.taper);
assert(max(abs(sum(abs(W).^2,1)-1))<1e-12);
tables_path=fullfile(parent,'all_data_fit','results','calibration_tables.mat');
tables=load(tables_path,'baseline_db','residual_db','commands');
az=-180:180;el=-90:90;
[AZ,EL]=meshgrid(az,el);
u=[cosd(EL(:)).'.*cosd(AZ(:)).';cosd(EL(:)).'.*sind(AZ(:)).';sind(EL(:)).'];
calibrated=qd_arrayant([]);calibrated.set_grid(az*pi/180,el*pi/180,0);
calibrated.no_elements=K;calibrated.center_frequency=M.fc;
calibrated.name='allData1970_interpolated_calibrated_beam_outputs';
Fa=zeros(numel(el),numel(az),K);Fb=Fa;power_grid=zeros(size(Fa));
corner_commands=cell(1,K);corner_weights=cell(1,K);
sign_j=-1;if strcmpi(polarization,'LHCP'),sign_j=1;end
for k=1:K
    scan=min(theta(k),50);p0=floor(scan/10)*10;p1=min(p0+10,50);
    if p0==p1,wp=0;else,wp=(scan-p0)/(p1-p0);end
    a0=floor(phi(k)/30)*30;a1=mod(a0+30,360);wa=(phi(k)-a0)/30;
    commands=[p0 a0;p0 a1;p1 a0;p1 a1];weights=[(1-wp)*(1-wa),(1-wp)*wa,wp*(1-wa),wp*wa];
    power=zeros(1,size(u,2));
    for corner=1:4
        if weights(corner)==0,continue;end
        source=dirvec(commands(corner,1),commands(corner,2));
        Q=align_rotation(source,U(:,k));query=Q.'*u;
        tq=acosd(max(-1,min(1,query(3,:))));aq=mod(atan2d(query(2,:),query(1,:)),360);
        tq(tq>90 & tq<90+1e-8)=90; % numerical boundary, not rear extrapolation
        idx=find(tables.commands(:,1)==commands(corner,1)&tables.commands(:,2)==commands(corner,2));
        assert(numel(idx)==1);
        gain=tables.baseline_db(:,:,idx)+tables.residual_db(:,:,idx);
        pg=10.^(gain/10);pg=[pg pg(:,1)];
        power=power+weights(corner)*interp2(0:360,(0:90).',pg,aq,tq,'linear',0);
    end
    power(u(3,:)<-1e-12)=0; % unmeasured physical rear hemisphere
    power=power*10^(gain_offset_dbi/10);
    % Retain the fitted 502-element geometric phase. This phase is assumed,
    % not measured; empirical calibration controls the combined-beam power.
    phase=zeros(1,size(u,2));
    for start=1:2048:size(u,2)
        ids=start:min(start+2047,size(u,2));
        af=M.taper(:).'*exp(2j*pi*M.fc/299792458*(M.p_site.'*(u(:,ids)-U(:,k))));
        phase(ids)=angle(af);phase(ids(abs(af)<1e-12))=0;
    end
    field=reshape(sqrt(max(power,0)).*exp(1j*phase),size(AZ));
    Fa(:,:,k)=field/sqrt(2);Fb(:,:,k)=sign_j*1j*Fa(:,:,k);
    power_grid(:,:,k)=reshape(power,size(AZ));
    corner_commands{k}=commands;corner_weights{k}=weights;
end
calibrated.Fa=Fa;calibrated.Fb=Fb;
calibrated.element_position=zeros(3,K);calibrated.coupling=eye(K);
info=struct('fc',M.fc,'N',N,'K',K,'p_canonical',M.p_site,'p_world',R*M.p_site, ...
    'rotation',R,'orientation',[pi;0;0],'Wq',W,'U_canonical',U,'U_world',U_world, ...
    'scan_theta_deg',theta,'scan_phi_deg',phi,'extrapolated',theta>50+1e-10, ...
    'corner_commands',{corner_commands},'corner_weights',{corner_weights}, ...
    'power_grid',power_grid,'az_grid_deg',az,'el_grid_deg',el, ...
    'polarization',upper(polarization),'gain_offset_dbi',gain_offset_dbi, ...
    'model_path',model_path,'calibration_path',tables_path);
end

function u=dirvec(theta,phi)
u=[sind(theta)*cosd(phi);sind(theta)*sind(phi);cosd(theta)];
end

function Q=align_rotation(source,target)
v=cross(source,target);c=max(-1,min(1,dot(source,target)));
assert(c>-1+1e-10,'Antipodal transport is undefined');
V=[0 -v(3) v(2);v(3) 0 -v(1);-v(2) v(1) 0];
Q=eye(3)+V+V*V/(1+c);
end
