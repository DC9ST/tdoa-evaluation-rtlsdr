function [ reliability ] = corr_reliability( correlation )
%corr_reliability calculates reliability 0..1 of the correlation result
%   Looks for 2nd highest outstanding PEAK (not sample!) and returns  "1 - amplitude percentage in comparison to main peak"
%	e.g. main peak 10, 2nd peak 3 -> output: 0.7

    [corr_max, idx] = max(correlation);

    correlation_temp = correlation;
    correlation_temp(idx) = 0;

    bin_right_old = corr_max;

    for ii=idx+1:1:length(correlation)          % delete all decreasing values on right side
        if (correlation_temp(ii) < bin_right_old)
            bin_right_old = correlation_temp(ii);
            correlation_temp(ii) = 0;
        else
            break;
        end
    end
    
    
    bin_left_old = corr_max;

    for ii=idx-1:-1:1                           % delete all decreasing values on left side
        if (correlation_temp(ii) < bin_left_old)
            bin_left_old = correlation_temp(ii);
            correlation_temp(ii) = 0;
        else
            break;
        end
    end

    peak2 = max(correlation_temp);              % look for 2nd (real) peak 
    reliability = 1-(peak2 / corr_max);
end


