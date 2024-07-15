% Parâmetros
EbN0 = 0:0.5:20;  % Valores de SNR em dB

% 64-QAM
[ser_64qam_theory, ber_64qam_theory] = ser_qam_theory(64, EbN0);
% 16-QAM
[ser_16qam_theory, ber_16qam_theory] = ser_qam_theory(16, EbN0);
% 8-PSK
[ser_8psk_theory, ber_8psk_theory] = ser_psk_theory(8, EbN0);
% QPSK
[ser_qpsk_theory, ber_qpsk_theory] = ser_psk_theory(4, EbN0);
% BPSK
EbN0_linear = 10.^(EbN0/10);
ser_bpsk_theory = qfunc(sqrt(2*EbN0_linear));
ber_bpsk_theory = ser_bpsk_theory;

% SER
figure;
semilogy(EbN0, ser_64qam_theory, 'o-', 'DisplayName', 'Teorico 64-QAM');
hold on;
semilogy(EbN0, ser_16qam_theory, 'o-', 'DisplayName', 'Teorico 16-QAM');
hold on;
semilogy(EbN0, ser_8psk_theory, 'o-', 'DisplayName', 'Teorico 8-PSK');
hold on;
semilogy(EbN0, ser_qpsk_theory, 'o-', 'DisplayName', 'Teorico QPSK');
hold on;
semilogy(EbN0, ser_bpsk_theory, 'o-', 'DisplayName', 'Teorico BPSK');
hold on;
ylim([10.^(-6) 10.^(-1)]);
xlabel('E_b/N_0 (dB)');
ylabel('Taxa de Erro de Símbolo (SER)');
legend show;
title('Taxa de Erro de Símbolo');
grid on;

% BER
figure;
semilogy(EbN0, ber_64qam_theory, 'o-', 'DisplayName', 'Teorico 64-QAM');
hold on;
semilogy(EbN0, ber_16qam_theory, 'o-', 'DisplayName', 'Teorico 16-QAM');
hold on;
semilogy(EbN0, ber_8psk_theory, 'o-', 'DisplayName', 'Teorico 8-PSK');
hold on;
semilogy(EbN0, ber_qpsk_theory, 'o-', 'DisplayName', 'Teorico QPSK');
hold on;
semilogy(EbN0, ber_bpsk_theory, 'o-', 'DisplayName', 'Teorico BPSK');
hold on;
ylim([10.^(-6) 10.^(0)]);
xlabel('E_b/N_0 (dB)');
ylabel('Taxa de Erro de Bit (BER)');
legend show;
title('Taxa de Erro de Bit');
grid on;



function [ser, ber] = ser_qam_theory(M, EbN0)
    EbN0_linear = 10.^(EbN0 / 10);
    ser = 4 * (1 - 1/sqrt(M)) * qfunc(sqrt(3 * (log2(M)/(M - 1)) * EbN0_linear));
    ber = ser / log2(M);
end

function [ser, ber] = ser_psk_theory(M, EbN0)
    EbN0_linear = 10.^(EbN0 / 10);
    ser = 2 * qfunc(sin(pi/M) * sqrt(2 * log2(M) * EbN0_linear));
    ber = ser / log2(M);
end
