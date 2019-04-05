function [ dist ] = dist_latlong( lat1, long1, lat2, long2, ref_lat, ref_long )
%DIST_LATLONG calculates distance between two specified points in meters
%accurate only in the (wide) area around the geodetic reference point ref_lat/_long

    [x1, y1] = latlong2xy( lat1, long1, ref_lat, ref_long );
    [x2, y2] = latlong2xy( lat2, long2, ref_lat, ref_long );

    dist = 1000 * sqrt( (x1-x2)^2 + (y1-y2)^2 );
end

