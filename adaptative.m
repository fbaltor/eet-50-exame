% Parâmetros
num_bits = 1000000; % Número de símbolos a serem transmitidos
data = randi([0 1], num_bits, 1);
EbN0 = 0:0.5:20;  % Valores de Eb/N0 em dB
Pe = [10.^(-2), 10.^(-3), 10.^(-4)];

% for i = 1:length(Pe)
%     plot_adaptative(data, EbN0, Pe(i), 'ber');
% end

plot_adaptative(data, EbN0, Pe(2), 'ser', [10.^(-7) 10.^(0)])
plot_adaptative(data, EbN0, Pe(2), 'ber', [10.^(-7) 10.^(0)])


function plot_adaptative(data, EbN0, pe, error_type, vlim)
    if nargin < 5
        vlim = [10.^(-6) 10.^(0)]; 
    end

    er = zeros(1, length(EbN0));
    
    for i = 1:length(EbN0)
        er(i) = metrics_adaptative(data, pe, EbN0(i), error_type);
    end
    
    er_64qam_theory = er_qam_theory(64, EbN0, error_type); % 64-QAM
    er_16qam_theory = er_qam_theory(16, EbN0, error_type); % 16-QAM
    er_8psk_theory = er_psk_theory(8, EbN0, error_type); % 8-PSK
    er_qpsk_theory = er_psk_theory(4, EbN0, error_type); % QPSK
    er_bpsk_theory = er_psk_theory(2, EbN0); % BPSK
    
    figure;
    semilogy(EbN0, er_64qam_theory, 'o-', 'DisplayName', 'Teorico 64-QAM');
    hold on;
    semilogy(EbN0, er_16qam_theory, 'o-', 'DisplayName', 'Teorico 16-QAM');
    hold on;
    semilogy(EbN0, er_8psk_theory, 'o-', 'DisplayName', 'Teorico 8-PSK');
    hold on;
    semilogy(EbN0, er_qpsk_theory, 'o-', 'DisplayName', 'Teorico QPSK');
    hold on;
    semilogy(EbN0, er_bpsk_theory, 'o-', 'DisplayName', 'Teorico BPSK');
    hold on;
    semilogy(EbN0, er, 'o-', 'DisplayName', 'Adaptativo');
    hold on;
    yline(pe, 'r', 'LineWidth', 1.2, 'DisplayName', 'Erro limite');
    ylim(vlim);
    xlabel('E_b/N_0 (dB)');

    if error_type == 'ser'
        text1 = "Taxa de Erro de Simbolo (SER)";
        text2 = "Simbolo";
    else
        text1 = "Taxa de Erro de Bit (BER)";
        text2 = "Bit";
    end

    ylabel(text1);
    legend show;
    title("Taxa de Erro de " + text2);
    grid on;
end

function er = metrics_adaptative(data_bits, pe, EbN0, error_type)
    if nargin < 4
        error_type = 'ser';
    end

    if er_qam_theory(64, EbN0, error_type) < pe
        er = metrics_qam(64, data_bits, EbN0, error_type);

    elseif er_qam_theory(16, EbN0, error_type) < pe
       er = metrics_qam(16, data_bits, EbN0, error_type);

    elseif er_psk_theory(8, EbN0, error_type) < pe
        er = metrics_psk(8, data_bits, EbN0, error_type);

    elseif er_psk_theory(4, EbN0, error_type) < pe
        er = metrics_psk(4, data_bits, EbN0, error_type);

    else
        er = metrics_psk(2, data_bits, EbN0, error_type);
    end
end

function er = metrics_psk(M, data_bits, EbN0, error_type)
    if nargin < 4
        error_type = 'ser';
    end

    k = log2(M);
    num_pad_bits = k - mod(length(data_bits), k);
    data_bits = [data_bits; zeros(num_pad_bits, 1)];
    data_symbols = bit2int(data_bits, k);
    tx_symbols = pskmod(data_symbols, M);
    snr = EbN0 + 10*log10(k);
    rx_symbols = awgn(tx_symbols, snr, 'measured');
    data_out_symbols = pskdemod(rx_symbols, M);
    data_out_bits = int2bit(data_out_symbols, k);
    [numErrors, ser] = symerr(data_symbols, data_out_symbols);
    [numErrors, ber] = biterr(data_bits, data_out_bits);

    if error_type == 'ser'
        er = ser;
    elseif error_type == 'ber'
        er = ber;
    end
end

function er = metrics_qam(M, data_bits, EbN0, error_type)
    if nargin < 4
        error_type = 'ser';
    end

    k = log2(M);
    num_pad_bits = k - mod(length(data_bits), k);
    data_bits = [data_bits; zeros(num_pad_bits, 1)];
    data_symbols = bit2int(data_bits, k);
    tx_symbols = qammod(data_symbols, M, 'UnitAveragePower', true);
    snr = EbN0 + 10*log10(k);
    rx_symbols = awgn(tx_symbols, snr, 'measured');
    data_out_symbols = qamdemod(rx_symbols, M, 'UnitAveragePower', true);
    data_out_bits = int2bit(data_out_symbols, k);
    [numErrors, ser] = symerr(data_symbols, data_out_symbols);
    [numErrors, ber] = biterr(data_bits, data_out_bits);

    if error_type == 'ser'
        er = ser;
    elseif error_type == 'ber'
        er = ber;
    end
end

function er = er_qam_theory(M, EbN0, error_type)
    if nargin < 3
        error_type = 'ser';
    end

    EbN0_linear = 10.^(EbN0 / 10);
    ser = 4 * (1 - 1/sqrt(M)) * qfunc(sqrt(3 * (log2(M)/(M - 1)) * EbN0_linear));
    ber = ser / log2(M);

    if error_type == 'ser'
        er = ser;
    elseif error_type == 'ber'
        er = ber;
    end
end

function er = er_psk_theory(M, EbN0, error_type)
    EbN0_linear = 10.^(EbN0 / 10);

    if M == 2
        er = qfunc(sqrt(2*EbN0_linear)); % BPSK
        return
    end

    if nargin < 3
        error_type = 'ser';
    end
    
    ser = 2 * qfunc(sin(pi/M) * sqrt(2 * log2(M) * EbN0_linear));
    ber = ser / log2(M);

    if error_type == 'ser'
        er = ser;
    elseif error_type == 'ber'
        er = ber;
    end
end