function AE = elem_gain_dB(az, el)
% ELEM_GAIN_DB  3GPP TR 38.901 element pattern [dBi], boresight along +x.
%   az, el : azimuth / elevation [deg] in the array frame (scalars or same-size arrays)
phi   = mod(az + 180, 360) - 180;     % wrap to [-180, 180)
theta = 90 - el;                      % zenith angle
AH = -min(12*(phi/65).^2, 30);
AV = -min(12*((theta-90)/65).^2, 30);
AE = 8 - min(-(AH + AV), 30);
end
