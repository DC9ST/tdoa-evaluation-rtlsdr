function [ x,y ] = latlong2xy_kl( lat, long )
% converts lat/long coordinates to cartesian x,y
% output in km

ref_lat = 49.4; % geodetic reference point near Kaiserslautern (for plane approximation)
ref_long = 7.7;
earth_radius = 40074; % km

y = (lat - ref_lat)/360 * earth_radius;
x = (long - ref_long)/360 * cos(ref_lat*pi/180) * earth_radius;

end

