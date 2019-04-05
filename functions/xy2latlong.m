function [ lat, long ] = xy2latlong( x, y, ref_lat, ref_long )
% converts cartesian coordinates x,y (in km) to lat/long
% ref_lat and ref_long is the geodetic reference point for plane approximation of the earth surface

earth_circumf = 40074;

lat = (y * 360 / earth_circumf) + ref_lat;
long = ( (x*360) / (earth_circumf * cos(ref_lat*pi/180)) ) + ref_long;

end

