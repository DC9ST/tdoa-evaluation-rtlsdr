function [ iq_corr ] = correlate_iq( iq1, iq2, corr_strategy, smoothing_factor )
%correlate_iq correlates two complex iq signals according to the specified strategy
% iq1: first complex IQ signal
% iq2: second complex IQ signal
% corr_strategy ={'abs', 'dphase'}
%   'abs': use absolute value of iq signals
%   'dphase': use differential phase of iq signals
% smoothing_factor

    switch corr_strategy
        
        case 'abs'
            abs1 = remove_mean(abs(iq1));
            abs2 = remove_mean(abs(iq2));

            abs_corr = xcorr(abs1, abs2);
            
            ref1 = max(xcorr(abs1, abs1));
            ref2 = max(xcorr(abs2, abs2));
            cor_max = max(abs_corr);
            disp(['abs cross-correlation: max (peak) ' num2str(cor_max,3) ', autocorr1 max ' num2str(ref1,3) ', autocorr2 max ' num2str(ref2,3) ', => ' num2str(100*2*cor_max / (ref1+ref2)) '%']);
            
            if smoothing_factor ~= 0 
                abs_corr = smooth(abs(abs_corr), smoothing_factor);
            end
            
            abs_corr1 = abs_corr ./ max(abs_corr); %normalize
            iq_corr = abs_corr1;
                        
        case 'dphase'
            d_phase1 = diff(unwrap(angle(iq1)));
            d_phase2 = diff(unwrap(angle(iq2)));
            
            d_phase1 = [ 0; d_phase1(1:length(d_phase1)) ];  % append a zero to match the length of the abs correlation (size fit)
            d_phase2 = [ 0; d_phase2(1:length(d_phase2)) ];

            d_phase1 = remove_mean(d_phase1);
            d_phase2 = remove_mean(d_phase2);        
            
            d_phase_corr = xcorr(d_phase1, d_phase2);

            ref1 = max(xcorr(d_phase1, d_phase1));
            ref2 = max(xcorr(d_phase2, d_phase2));
            cor_max = max(d_phase_corr);
            disp(['dphase cross-correlation, max (peak) ' num2str(cor_max,3) ', autocorr1 max ' num2str(ref1,3) ', autocorr2 max ' num2str(ref2,3) ', => ' num2str(100*2*cor_max / (ref1+ref2)) '%'  ]);
            
            if smoothing_factor ~= 0 
                d_phase_corr = smooth(abs(d_phase_corr), smoothing_factor);
            end
            
            d_phase_corr1 = d_phase_corr ./ max(d_phase_corr); %normalize
            iq_corr = d_phase_corr1;
            
        otherwise
            error('correlate_iq.m: correlation strategy not supported (choose abs or dphase)');
    end 

end

