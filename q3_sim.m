% Parâmetros
numSymbols = 1000000; % Número de símbolos a serem transmitidos

data = randi([0 1], numSymbols, 1);

snr = 0:0.5:20;  % Valores de SNR em dB
snr_linear = 10.^(snr/10);

% 64-QAM
ser_64qam = zeros(1, length(snr));
ser_64qam_theory = 4 * (1 - 1/sqrt(64)) * qfunc(sqrt((3 * log2(64)/(M - 1)) * snr_linear));

% 16-QAM
ser_16qam = zeros(1, length(snr));
ser_16qam_theory = 2 * (1 - 1/16) * qfunc(sqrt(6 * (log2(16)/(M^2 - 1)) * snr_linear));

% 8-PSK
ser_8psk = zeros(1, length(snr));
ser_8psk_theory = 2 * qfunc(sqrt(2 * log2(8) * ))

ser_qpsk = zeros(1, length(snr));  % SER simulado
ser_bpsk = zeros(1, length(snr));  % SER simulado


for i = 1:length(snr)

    % 64-QAM
    modData_64qam = qammod(data, 64, 'UnitAveragePower', true);
    receivedData_64qam = awgn(modData_64qam, snr(i), 'measured');
    demodData_64qam = qamdemod(receivedData_64qam, 64, 'UnitAveragePower', true);
    [numErrors, ratio] = symerr(data, demodData_64qam);
    ser_64qam(i) = numErrors / numSymbols;


    % 16-QAM
    modData_16qam = qammod(data, 16, 'UnitAveragePower', true);
    receivedData_16qam = awgn(modData_16qam, snr(i), 'measured');
    demodData_16qam = qamdemod(receivedData_16qam, 16, 'UnitAveragePower', true);
    [numErrors, ratio] = symerr(data, demodData_16qam);
    ser_16qam(i) = numErrors / numSymbols;

    % 8-PSK
    modData_8psk = pskmod(data, 8);
    receivedData_8psk = awgn(modData_8psk, snr(i), 'measured');
    demodData_8psk = pskdemod(receivedData_8psk, 8);
    [numErrors, ratio] = symerr(data, demodData_8psk);
    ser_8psk(i) = numErrors / numSymbols;

    % QPSK
    modData_qpsk = pskmod(data, 4);
    receivedData_qpsk = awgn(modData_qpsk, snr(i), 'measured');
    demodData_qpsk = pskdemod(receivedData_qpsk, 4);
    [numErrors, ratio] = symerr(data, demodData_qpsk);
    ser_qpsk(i) = numErrors / numSymbols;

    % BPSK
    modData_bpsk = pskmod(data, 2);
    receiveData_bpsk = awgn(modData_bpsk, snr(i), 'measured');
    demodData_bpsk = pskdemod(receiveData_bpsk, 2);
    [numErrors, ratio] = symerr(data, demodData_bpsk);
    ser_bpsk(i) = numErrors / numSymbols;
end

figure;
semilogy(snr, ser_64qam, 'o-', 'DisplayName', 'Simulação 64-QAM');
hold on;
semilogy(snr, ser_16qam, 'o-', 'DisplayName', 'Simulação 16-QAM');
hold on;
semilogy(snr, ser_8psk, 'o-', 'DisplayName', 'Simulação 8-PSK');
hold on;
semilogy(snr, ser_qpsk, 'o-', 'DisplayName', 'Simulação QPSK');
hold on;
semilogy(snr, ser_bpsk, 'o-', 'DisplayName', 'Simulação BPSK');
hold on;
xlabel('E_b/N_0 (dB)');
ylabel('Taxa de Erro de Símbolo (SER)');
legend show;
title('Taxa de Erro de Símbolo');
grid on;