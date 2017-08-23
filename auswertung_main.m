% =========================================================================
%  Experimental Evaluation Script for RTL-SDR based TDOA
%  DC9ST, 2017
% =========================================================================

clear;
clc;
close all;

% adds subfolder with functions to PATH
[p,n,e] = fileparts(mfilename('fullpath'));
addpath([p '\functions'])


%% Specify Parameters

% RX and Ref TX Position
rx1_lat = 49.441; % RX 1
rx1_long = 7.767;

rx2_lat = 49.422; % RX 2
rx2_long = 7.739;

rx3_lat = 49.425; % RX 3
rx3_long = 7.756;

tx_ref_lat = 49.45962; % Referenz: Rotenberg DAB
tx_ref_long = 7.77116;

% signal processing parameters
signal_bandwidth_khz = 0;  % 400, 200, 40, 12, 0(no)
smoothing_factor = 0;
corr_type = 'dphase';  %'abs' or 'dphase'
interpol_factor = 0;

% 1: show correlation plots
% 2: show also input spcetrograms and spectra of input meas
% 3: show also before and after filtering
report_level = 1;

% heatmap
heatmap_resolution = 400; % resolution for heatmap points
heatmap_threshold = 0.1;  % heatmap point with lower mag are suppressed for html output

% IQ Data Files
file_identifier = 'test.dat';

dateiname1 = ['recorded_data\1_' file_identifier];
dateiname2 = ['recorded_data\2_' file_identifier];
dateiname3 = ['recorded_data\3_' file_identifier];

% known signal path differences between two RXes to Ref (sign of result is important!)
rx_distance_diff12 = dist_latlong_kl(tx_ref_lat, tx_ref_long, rx1_lat, rx1_long) - dist_latlong_kl(tx_ref_lat, tx_ref_long, rx2_lat, rx2_long); % (Ref to RX1 - Ref to RX2) in meters
rx_distance_diff13 = dist_latlong_kl(tx_ref_lat, tx_ref_long, rx1_lat, rx1_long) - dist_latlong_kl(tx_ref_lat, tx_ref_long, rx3_lat, rx3_long); % (Ref to RX1 - Ref to RX3) in meters
rx_distance_diff23 = dist_latlong_kl(tx_ref_lat, tx_ref_long, rx2_lat, rx2_long) - dist_latlong_kl(tx_ref_lat, tx_ref_long, rx3_lat, rx3_long); % (Ref to RX2 - Ref to RX3) in meters

% distance between two RXes in meters
rx_distance12 = dist_latlong_kl(rx1_lat, rx1_long, rx2_lat, rx2_long);
rx_distance13 = dist_latlong_kl(rx1_lat, rx1_long, rx3_lat, rx3_long);
rx_distance23 = dist_latlong_kl(rx2_lat, rx2_long, rx3_lat, rx3_long);

%% Read Signals from File
disp('______________________________________________________________________________________________');
disp('READ DATA FROM FILES');
signal1 = read_file_iq(dateiname1);
signal2 = read_file_iq(dateiname2);
signal3 = read_file_iq(dateiname3);

if (report_level > 1)
    % display raw signals
    num_samples_total = length(signal1);
    inphase1 = real(signal1);
    quadrature1 = imag(signal1);
    inphase2 = real(signal2);
    quadrature2 = imag(signal2);
    inphase3 = real(signal3);
    quadrature3 = imag(signal3);
    
    figure;
    subplot(3,1,1);
    plot(1:num_samples_total, inphase1(1:num_samples_total), 1:num_samples_total, quadrature1(1:num_samples_total));
    title('raw RX 1: I and Q');
    subplot(3,1,2);
    plot(1:num_samples_total, inphase2(1:num_samples_total), 1:num_samples_total, quadrature2(1:num_samples_total));
    title('raw RX 2: I and Q');    
    subplot(3,1,3);
    plot(1:num_samples_total, inphase3(1:num_samples_total), 1:num_samples_total, quadrature3(1:num_samples_total));
    title('raw RX 3: I and Q');
end

if (report_level > 1)
    % calculate and show spectrogram
    nfft = 256;
    overlap = -1000;
    
    figure;
    subplot(4,2,1);
    complex_signal = detrend(signal1);
    [S,F,T,P] = spectrogram(complex_signal, nfft, overlap, nfft, 2e6 );
    spectrum = fftshift(fliplr(10*log10(abs(P))'), 2);
    for i=1:nfft
        spectrum(:,i) = smooth(spectrum(:,i),9);
    end
    surf(T,F, spectrum', 'edgecolor', 'none');
    axis tight;
    view(0,90);
    title('RX 1');
    xlabel('time');
    ylabel('frequency');

    subplot(4,2,3);
    complex_signal = detrend(signal2);
    [S,F,T,P] = spectrogram(complex_signal, nfft, overlap, nfft, 2e6 );
    spectrum = fftshift(fliplr(10*log10(abs(P))'), 2);
    for i=1:nfft
        spectrum(:,i) = smooth(spectrum(:,i),9);
    end
    surf(T,F, spectrum', 'edgecolor', 'none');
    axis tight;
    view(0,90);
    title('RX 2');
    xlabel('time');
    ylabel('frequency');

    
    subplot(4,2,5);
    complex_signal = detrend(signal3);
    [S,F,T,P] = spectrogram(complex_signal, nfft, overlap, nfft, 2e6 );
    spectrum = fftshift(fliplr(10*log10(abs(P))'), 2);
    for i=1:nfft
        spectrum(:,i) = smooth(spectrum(:,i),9);
    end
    surf(T,F, spectrum', 'edgecolor', 'none');
    axis tight;
    view(0,90);
    title('RX 3');
    xlabel('time');
    ylabel('frequency');
    
    % display spectrum
    spectrum_smooth_factor  = 201; 
    subplot(4,2,2);
    spectrum_single1 = 10*log10(abs(fftshift(fft(signal1(1.7e6 : 1.7e6 + 2^18)))));
    spectrum_single1 = smooth(spectrum_single1, spectrum_smooth_factor);
    plot(spectrum_single1);
    title('Measurement RX 1');
    grid;

    subplot(4,2,4);
    spectrum_single2 = 10*log10(abs(fftshift(fft(signal2(1.7e6 : 1.7e6 + 2^18)))));
    spectrum_single2 = smooth(spectrum_single2, spectrum_smooth_factor);
    plot(spectrum_single2);
    title('Measurement RX 2');
    grid;

    subplot(4,2,6);
    spectrum_single3 = 10*log10(abs(fftshift(fft(signal3(1.7e6 : 1.7e6 + 2^18)))));
    spectrum_single3 = smooth(spectrum_single3, spectrum_smooth_factor);
    plot(spectrum_single3);
    title('Measurement RX 3');
    grid;

    subplot(4,2,7:8);
    freq_axis = -(length(spectrum_single1)/2) : 1 : ((length(spectrum_single1)/2)-1);
    plot(freq_axis, spectrum_single1, freq_axis, spectrum_single2, freq_axis, spectrum_single3);
    title('Measurement Signal RX 1,2 & 3');
    grid;
end;



%% Calculate TDOA
disp(' ');
disp('______________________________________________________________________________________________');
disp('CORRELATION 1 & 2');
[doa_meters12, doa_samples12, reliability12 ] = tdoa2(signal1, signal2, rx_distance_diff12, rx_distance12, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor);

disp(' ');
disp('______________________________________________________________________________________________');
disp('CORRELATION 1 & 3');
[doa_meters13, doa_samples13, reliability13 ] = tdoa2(signal1, signal3, rx_distance_diff13, rx_distance13, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor);

disp(' ');
disp('______________________________________________________________________________________________');
disp('CORRELATION 2 & 3');
[doa_meters23, doa_samples23, reliability23 ] = tdoa2(signal2, signal3, rx_distance_diff23, rx_distance23, smoothing_factor, corr_type, report_level, signal_bandwidth_khz, interpol_factor);


%% Generate html map
disp(' ');
disp('______________________________________________________________________________________________');
disp('GENERATE HYPERBOLAS');

[points_lat1, points_long1] = gen_hyperbola(doa_meters12, rx1_lat, rx1_long, rx2_lat, rx2_long);
[points_lat2, points_long2] = gen_hyperbola(doa_meters13, rx1_lat, rx1_long, rx3_lat, rx3_long);
[points_lat3, points_long3] = gen_hyperbola(doa_meters23, rx2_lat, rx2_long, rx3_lat, rx3_long);

disp(' ');
disp('______________________________________________________________________________________________');
disp('GENERATE HTML');
rx_lat_positions  = [rx1_lat   rx2_lat   rx3_lat ];
rx_long_positions = [rx1_long  rx2_long  rx3_long];

hyperbola_lat_cell  = {points_lat1,  points_lat2, points_lat3};
hyperbola_long_cell = {points_long1, points_long2, points_long3};

[heatmap_long, heatmap_lat, heatmap_mag] = create_heatmap_kl(doa_meters12, doa_meters13, doa_meters23, rx1_lat, rx1_long, rx2_lat, rx2_long, rx3_lat, rx3_long, heatmap_resolution); % generate heatmap
heatmap_cell = {heatmap_long, heatmap_lat, heatmap_mag};

create_html_file( ['ergebnisse/map_' file_identifier '_' corr_type '_interp' num2str(interpol_factor) '_bw' int2str(signal_bandwidth_khz) '_smooth' int2str(smoothing_factor) '.html'], rx_lat_positions, rx_long_positions, hyperbola_lat_cell, hyperbola_long_cell, heatmap_cell, heatmap_threshold);

disp('______________________________________________________________________________________________');