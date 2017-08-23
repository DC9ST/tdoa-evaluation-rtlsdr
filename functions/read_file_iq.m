function [ iq_signal ] = read_file_iq( filename )
%read_file_tdoa reads a file with IQ data (e.g. from rtl_sdr)

    disp('read_file_iq');

    % daten vom 1. RX
    disp(['IQ read from data file = ' filename]);
    fileID = fopen(filename);
    a = fread(fileID);
    fclose(fileID);

    inphase1 = a(1:2:end) -128;
    quadrature1 = a(2:2:end) -128;
    disp(['successfully read ' int2str(length(inphase1)) ' samples']);

    % complex representation
    iq_signal = inphase1 + 1i.*quadrature1;

end

