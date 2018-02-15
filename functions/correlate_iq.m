function [ iq_corr ] = correlate_iq( iq1, iq2, corr_strategy, smoothing_factor )
%correlate_iq correlates two complex iq signal according to the specified strategy
% iq1: first complex IQ signal
% iq2: second complex IQ signal
% corr_strategy ={'abs', 'dphase'}
%   'abs': use absolute value of iq signals (with detrend)
%   'dphase': use differential phase of iq signals (with detrend and one zero inserted at beginning for size fit)
% smoothing_factor

    switch corr_strategy
        
        case 'abs'
            abs1_detrend = detrend(abs(iq1));
            abs2_detrend = detrend(abs(iq2));

            abs_detrend_corr = xcorr(abs1_detrend, abs2_detrend);
            
            ref1 = max(xcorr(abs1_detrend, abs1_detrend));
            ref2 = max(xcorr(abs2_detrend, abs2_detrend));
            cor_max = max(abs_detrend_corr);
            disp(['abs Korrelation, max ' num2str(cor_max) ', ref1 = ' num2str(ref1) ', ref2= ' num2str(ref2) ', => ' num2str(100*2*cor_max / (ref1+ref2)) '%']);
            
            if smoothing_factor ~= 0 
                abs_detrend_corr = smooth(abs(abs_detrend_corr), smoothing_factor);
            end
            
            abs_detrend_corr1 = abs_detrend_corr ./ max(abs_detrend_corr); %normalize
            iq_corr = abs_detrend_corr1;
                        
        case 'dphase'
            d_phase1 = diff(unwrap(angle(iq1)));
            d_phase2 = diff(unwrap(angle(iq2)));
            
            d_phase1 = [ 0; d_phase1(1:length(d_phase1)) ];  % append a zero to match the length of the abs correlation
            d_phase2 = [ 0; d_phase2(1:length(d_phase2)) ];

            d_phase1_detrend = detrend(d_phase1);
            d_phase2_detrend = detrend(d_phase2);        
            
            d_phase_detrend_corr = xcorr(d_phase1_detrend, d_phase2_detrend);

            ref1 = max(xcorr(d_phase1_detrend,d_phase1_detrend));
            ref2 = max(xcorr(d_phase2_detrend,d_phase2_detrend));
            cor_max = max(d_phase_detrend_corr);
            disp(['dphase Korrelation, max ' num2str(cor_max) ', ref1 = ' num2str(ref1) ', ref2= ' num2str(ref2) ', => ' num2str(100*2*cor_max / (ref1+ref2)) '%'  ]);
            
            if smoothing_factor ~= 0 
                d_phase_detrend_corr = smooth(abs(d_phase_detrend_corr), smoothing_factor);
            end
            
            d_phase_detrend_corr1 = d_phase_detrend_corr ./ max(d_phase_detrend_corr); %normalize
            iq_corr = d_phase_detrend_corr1;
            
        otherwise
            error('correlate_iq.m: correlation strategy not supported (choose abs or dphase)');
    end 

end

