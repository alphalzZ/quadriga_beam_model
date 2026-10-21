function [W, U] = build_beam_codebook(tilt_deg, azi_deg, p, fc)
% BUILD_BEAM_CODEBOOK  Receive beam weights for a planar array (GNU Octave / MATLAB)
%   tilt_deg : downtilt angle [deg], positive = below horizon
%   azi_deg  : azimuth [deg] in the array frame (0 = broadside, +x)
%   p        : 3 x N element positions [m]
%   fc       : carrier frequency [Hz]
%   W        : N x K, unit-norm weight vector per beam, w = exp(+j k0 p.'u)/sqrt(N)
%   U        : 3 x K, unit direction vectors of the beams
c0 = 299792458;
k0 = 2*pi*fc/c0;
ele = -tilt_deg(:).';                 % elevation above horizon
az  = azi_deg(:).';
U   = [cosd(ele).*cosd(az); cosd(ele).*sind(az); sind(ele)];
N   = size(p, 2);
W   = exp(1j*k0*(p.'*U)) / sqrt(N);
end
