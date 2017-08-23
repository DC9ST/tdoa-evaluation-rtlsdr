function [ heat_long, heat_lat, mse_doa ] = create_heatmap_kl( doa_meters12, doa_meters13, doa_meters23, rx1_lat, rx1_long, rx2_lat, rx2_long, rx3_lat, rx3_long, resolution )
%create_heatmap_kl Creates a heatmap for Kaiserslautern based on mean squared error

%    returns: 
%    heat_long: longitudes of heatmap points
%    heat_lat: latitudes of heatmap points
%    mse_doa: heatmap magnitudes

    disp('creating heatmap... ');

    num_points = resolution; % points in one dimension (creates squared area)

    start_lat  = 49.40;
    stop_lat   = 49.50;
    %stop_lat   = 49.55;
    %step_lat   = (stop_lat - start_lat) / num_points
    
    start_long = 7.65;
    stop_long  = 7.87;
    %step_long  = (stop_long - start_long) / num_points
    
    heat_lat  = linspace(start_lat,  stop_lat,  num_points);
    heat_long = linspace(start_long, stop_long, num_points);
    mse_doa = zeros(num_points, num_points);
    
    for lat_idx = 1:num_points
        for long_idx = 1:num_points
            % calculate mean squared error of current point in terms of tdoa

            % distance current point to receivers
            dist_to_rx1 = dist_latlong_kl( heat_lat(lat_idx), heat_long(long_idx), rx1_lat, rx1_long );
            dist_to_rx2 = dist_latlong_kl( heat_lat(lat_idx), heat_long(long_idx), rx2_lat, rx2_long );
            dist_to_rx3 = dist_latlong_kl( heat_lat(lat_idx), heat_long(long_idx), rx3_lat, rx3_long );
            
            % current doa in meters
            current_doa12 = dist_to_rx1 - dist_to_rx2;
            current_doa13 = dist_to_rx1 - dist_to_rx3;
            current_doa23 = dist_to_rx2 - dist_to_rx3;
            
            % error doa
            doa_error = (current_doa12 - doa_meters12)^2 + (current_doa13 - doa_meters13)^2 + (current_doa23 - doa_meters23)^2;
            mse_doa(long_idx, lat_idx) = doa_error;
        end
    end
    
    mse_doa = 1./mse_doa;
    
    mse_doa = mse_doa .* (1/max(max(mse_doa)));
    
    %mse_doa = ones(num_points, num_points);
    
    disp('creating heatmap done! ');
end

