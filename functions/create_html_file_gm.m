function [ output_args ] = create_html_file_gm( filename, rx_lat, rx_long, hyperbola_lat, hyperbola_long, heatmap_cell, heatmap_threshold )
%CREATE_HTML_FILE   Generates the html code for google maps showing RX
%                   positions and the hyperbolas with a heatmap
%
%   rx_lat:         list of latitudes for RX positions (for displaying receiver positions)
%   rx_long:        list of longitudes for RX positions
%   hyperbola_lat:  cell array of latitudes for RX positions (one cell array elements correspond to one hyperbola)
%   hyperbola_lat:  cell array of longitudes for RX positions (one cell array elements correspond to one hyperbola)
%   heatmap_cell:   cell array, that contains {1} = long vector, {2} = lat vector and {3} = magnitude of heatmap points vector

    disp('writing html for google maps... ');

    % consistency checks
    if ~iscell(heatmap_cell)
        error('Parameter heatmap_cell needs to be a cell array with long, lat and magnitude');
    end
    
    if size(heatmap_cell) ~= 3
        error('Parameter heatmap_cell needs to be a cell array with long, lat and magnitude');
    end

    
    if ~iscell(hyperbola_lat) || ~iscell(hyperbola_long)
        error('Parameter hyperbola_lat, hyperbola_long need to be cell arrays');
    end
    
    if length(hyperbola_lat) ~= length(hyperbola_long)
        error('Length of hyperbola latitude and longitude values do not match');
    end

    for ii=1:length(hyperbola_lat)
        if length(cell2mat(hyperbola_lat(ii))) ~= length(cell2mat(hyperbola_long(ii)))
            error(['Length of hyperbola latitude and longitude in cell array doesn not for hyperbola ' num2str(ii)]);
        end
    end 

    if size(rx_lat) ~= size(rx_long)
        error('Dimensions of rx latitude and longitude values do not match');
    end
    
    
    num_rx_positions = length(rx_lat);

    num_hyperbolas = length(hyperbola_lat);  % length of cell array

    num_hyperb_points = zeros(1, num_hyperbolas);
    for ii = 1:num_hyperbolas
        num_hyperb_points(ii) = length(hyperbola_lat{ii}); % length of each vector in the cell array    
    end
    

    %% generate html code
    disp(['write html to file: ' filename]);

    fid = fopen(filename,'w');
  
    fprintf(fid, [...
        '<!DOCTYPE html>\n'...
        '<html>\n'...
        '\n'...
        '<body>\n'...
        '\n'...
        '<div id="map" style="width:100%%;height:800px"></div>\n'...
        '\n'...
        '<script>\n'...
        'function myMap() { \n'...
        '  var myCenter = new google.maps.LatLng(' num2str(mean(rx_lat), 8) ', ' num2str(mean(rx_long), 8) '); \n'...
        '  var mapCanvas = document.getElementById("map"); \n'...
        '  var mapProp = {center: myCenter, zoom: 13, scaleControl: true}; \n'...
        '  var map = new google.maps.Map(mapCanvas, mapProp); \n'...
        ]);

    %% write RX positions to html file
    for ii = 1:num_rx_positions
        fprintf(fid, '\n');
        fprintf(fid, ['  var rx' num2str(ii) '_MarkerCenter = new google.maps.LatLng(' num2str(rx_lat(ii), 8) ', ' num2str(rx_long(ii), 8) '); \n']);
        fprintf(fid, ['  var marker_rx' num2str(ii) ' = new google.maps.Marker({position: rx' num2str(ii) '_MarkerCenter}); \n']);
        fprintf(fid, ['  marker_rx' num2str(ii) '.setMap(map); \n\n']);
    end;
    
    %% write hyperbolas
    %define hyperbola arrays, one for each hyperbola
    for ii = 1:num_hyperbolas
        fprintf(fid, ['  const hyperbola_' num2str(ii) ' = Array(' num2str(num_hyperb_points(ii)) ') \n']);    
    end;
    fprintf(fid, '\n');
    
    % fill arrays with data points
    for ii_hyperb = 1:num_hyperbolas
        hyperbola_lat_vector  = cell2mat(hyperbola_lat(ii_hyperb));
        hyperbola_long_vector = cell2mat(hyperbola_long(ii_hyperb));
        
        for ii_point = 1:num_hyperb_points(ii_hyperb)
           fprintf(fid, ['    hyperbola_' num2str(ii_hyperb) '[' num2str(ii_point) '] = new google.maps.LatLng(' num2str(hyperbola_lat_vector(ii_point),8) ',' num2str(hyperbola_long_vector(ii_point),8) '); \n']);
        end
        fprintf(fid, '\n');         
    end
    fprintf(fid, '\n');
    
    % collect points of the array for each hyperbola
    for ii_hyperb = 1:num_hyperbolas
        fprintf(fid, ['  var path_hyperbola_' num2str(ii_hyperb) ' = [']);
        
        for ii_point = 1:num_hyperb_points(ii_hyperb)
            fprintf(fid, ['hyperbola_' num2str(ii_hyperb) '[' num2str(ii_point) ']' ]);
            if ii_point < num_hyperb_points(ii_hyperb)
                 fprintf(fid, ', ');
            %else
            %    fprintf(fid,'];\n');
            end
        end
        fprintf(fid,'];\n');
    end
    
    fprintf(fid, '\n');
    
    for ii_hyperb = 1:num_hyperbolas
        fprintf(fid, ['  var flightPath' num2str(ii_hyperb) '= new google.maps.Polyline({\n'...
                      '    path: path_hyperbola_' num2str(ii_hyperb) ',\n' ...
                      '    strokeColor:"#FF0000",\n'...
                      '    strokeOpacity:0.6,\n'...
                      '    strokeWeight:3\n'...
                      '  });\n'...
                      '  \n'...
                      ]);
    end
    
    for ii_hyperb = 1:num_hyperbolas
        fprintf(fid, ['  flightPath' num2str(ii_hyperb) '.setMap(map); \n']);
    end
    
    
    %% add heatmap
    
    %  extract data
    heat_long  = heatmap_cell{1};
    heat_lat  = heatmap_cell{2};
    heat_mag  = heatmap_cell{3};
    
    if ((length(heat_long) ~= length(heat_lat)) || (length(heat_long) ~= length(heat_mag)))
        error('create_html_file.m: Length of heatmap vectors in cell array do not match');
    end
    
    heat_num_points = length(heat_long);
    
    fprintf(fid, '\n  var heatmapData = [\n');
            
    for ii_lat = 1:heat_num_points
        for ii_long = 1:heat_num_points
            if (heat_mag(ii_long, ii_lat) > heatmap_threshold) || ((ii_lat == heat_num_points) && (ii_long == heat_num_points))
                fprintf(fid, ['    {location: new google.maps.LatLng(' num2str(heat_lat(ii_lat)) ', ' num2str(heat_long(ii_long)) '), weight: ' num2str(heat_mag(ii_long, ii_lat)) '}']) ;
       
                if (ii_lat == heat_num_points) && (ii_long == heat_num_points)
                    fprintf(fid, '\n  ];\n'); 
                else
                    fprintf(fid, ',\n');
                end
            end
            
        end    
    end
    
    
    fprintf(fid, [...
        '  \n'...
        '  var heatmap = new google.maps.visualization.HeatmapLayer({\n'...
        '    data: heatmapData,\n'...
        '    radius: 10\n'...
	    '    });\n'...
        '\n'...
	    '  heatmap.setMap(map);\n'...
        ]);    
    
   
    
    
    
    %% footer
    
    fprintf(fid, [...
        '\n'...
        '} \n'...
        '\n\n'...
        '</script>\n'...
        '\n'...
        '<script src="https://maps.googleapis.com/maps/api/js?callback=myMap&libraries=visualization"></script>\n'...
        '\n'...
        '</body>\n'...
        '</html>\n'...
        '\n'...
        ]);

    fclose(fid);
    
    disp('writing html, done!');
end

