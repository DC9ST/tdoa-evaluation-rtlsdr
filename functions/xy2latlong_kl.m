function [ lat, long ] = xy2latlong_kl( x, y )
% converts cartesian coordinates x,y to lat/long

ref_lat = 49.4; % geodetic reference point near Kaiserslautern (for plane approximation)
ref_long = 7.7;
earth_circumf = 40074;

lat = (y * 360 / earth_circumf) + ref_lat;
long = ( (x*360) / (earth_circumf * cos(ref_lat*pi/180)) ) + ref_long;

end

