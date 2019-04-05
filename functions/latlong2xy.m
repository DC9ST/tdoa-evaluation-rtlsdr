function [ x,y ] = latlong2xy( lat, long, ref_lat, ref_long )
% converts lat/long coordinates to cartesian x,y, output in km
% ref_lat and ref_long is the geodetic reference point for plane approximation of the earth surface

earth_circumf = 40074; % km

y = (lat - ref_lat)/360 * earth_circumf;
x = (long - ref_long)/360 * cos(ref_lat*pi/180) * earth_circumf;

end

