function d = angle_diff_deg(a, b)
% ANGLE_DIFF_DEG  Smallest signed difference a - b in degrees, wrapped to [-180, 180).
d = mod(a - b + 180, 360) - 180;
end
