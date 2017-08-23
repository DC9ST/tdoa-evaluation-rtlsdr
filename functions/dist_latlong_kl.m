function [ dist ] = dist_latlong_kl( lat1, long1, lat2, long2 )
%DIST_LATLONG_KL calculates distance between two specified points in meters
%around Kaiserslautern

    [x1, y1] = latlong2xy_kl( lat1, long1 );
    [x2, y2] = latlong2xy_kl( lat2, long2 );

    dist = 1000 * sqrt( (x1-x2)^2 + (y1-y2)^2 );
end

