function [ signal_iq_filtered ] = filter_iq( signal_iq, signal_bandwidth_khz )
%filter_iq filters an iq signal with a FIR filter of the specified
%bandwidth, if bandwidth is not supported, no filtering is performed
% signal_bandwidth_khz: full bandwidth of signal (filter uses bandwidth/2),
% possible values: 400, 200 (UKW), 12 (12.5, NBFM channel), 40 

    switch signal_bandwidth_khz
        
        case 400        
            % Equiripple Lowpass filter designed using the FIRPM function.
            % f_pass = 200 kHz, f_stop = 300 kHz, Apass = 0.1, Astop = 60

            % All frequency values are in kHz.
            Fs = 2000;  % Sampling Frequency

            Fpass = 200;              % Passband Frequency
            Fstop = 300;              % Stopband Frequency
            Dpass = 0.0057563991496;  % Passband Ripple
            Dstop = 0.001;            % Stopband Attenuation
            dens  = 20;               % Density Factor

            % Calculate the order from the parameters using FIRPMORD.
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

            % Calculate the coefficients using the FIRPM function.
            b  = firpm(N, Fo, Ao, W, {dens});
            Hd = dfilt.dffir(b);
            
            %filter
            signal_iq_filtered = filter(Hd, signal_iq);
            disp('Signal filtered to 400 kHz');
            
        case 200
            % Equiripple Lowpass filter designed using the FIRPM function.
            % f_pass = 100 kHz, f_stop = 150 kHz, Apass = 0.5, Astop = 60

            % All frequency values are in kHz.
            Fs = 2000;  % Sampling Frequency

            Fpass = 100;             % Passband Frequency
            Fstop = 150;             % Stopband Frequency
            Dpass = 0.028774368332;  % Passband Ripple
            Dstop = 0.001;           % Stopband Attenuation
            dens  = 20;              % Density Factor

            % Calculate the order from the parameters using FIRPMORD.
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

            % Calculate the coefficients using the FIRPM function.
            b  = firpm(N, Fo, Ao, W, {dens});
            Hd = dfilt.dffir(b);
        
            %filter
            signal_iq_filtered = filter(Hd, signal_iq);
            disp('Signal filtered to 200 kHz');

        case 40     
            % Equiripple Lowpass filter designed using the FIRPM function.
            % f_pass = 20 kHz, f_stop = 100 kHz, Apass = 0.1, Astop = 60

            % All frequency values are in kHz.
            Fs = 2000;  % Sampling Frequency

            Fpass = 20;               % Passband Frequency
            Fstop = 100;              % Stopband Frequency
            Dpass = 0.0057563991496;  % Passband Ripple
            Dstop = 0.001;            % Stopband Attenuation
            dens  = 20;               % Density Factor

            % Calculate the order from the parameters using FIRPMORD.
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

            % Calculate the coefficients using the FIRPM function.
            b  = firpm(N, Fo, Ao, W, {dens});
            Hd = dfilt.dffir(b);
            
            %filter
            signal_iq_filtered = filter(Hd, signal_iq);
            disp('Signal filtered to 40 kHz');

        case 12
            % Equiripple Lowpass filter designed using the FIRPM function.
            % f_pass = 6.25 kHz, f_stop = 50 kHz, Apass = 0.5, Astop = 60
            
            % All frequency values are in kHz.
            Fs = 2000;  % Sampling Frequency

            Fpass = 6.25;            % Passband Frequency
            Fstop = 50;              % Stopband Frequency
            Dpass = 0.028774368332;  % Passband Ripple
            Dstop = 0.001;           % Stopband Attenuation
            dens  = 20;              % Density Factor

            % Calculate the order from the parameters using FIRPMORD.
            [N, Fo, Ao, W] = firpmord([Fpass, Fstop]/(Fs/2), [1 0], [Dpass, Dstop]);

            % Calculate the coefficients using the FIRPM function.
            b  = firpm(N, Fo, Ao, W, {dens});
            Hd = dfilt.dffir(b);
            
            %filter
            signal_iq_filtered = filter(Hd, signal_iq);
            disp('Signal filtered to 12.5 kHz');
            
        case 0
            signal_iq_filtered = signal_iq;
            disp('Signal not filtered');
            
        otherwise
            signal_iq_filtered = signal_iq;
            disp('Specified filter bandwidth invalid!!');
            disp('no filtering performed!');
    end    

end

