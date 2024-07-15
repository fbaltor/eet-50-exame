% Simulação de 4-PAM

% Parâmetros
N = 1e5;  % Número de símbolos
M = 4;  % Ordem da modulação
k = log2(M);  % Bits por símbolo
EbN0_dB = 0:10;  % Valores de Eb/N0 em dB
Pe_sim = zeros(1,length(EbN0_dB));  % SER simulado
simbolos=[-9 -3 3 9];

% Simulação Monte Carlo
for i = 1:length(EbN0_dB)
    EbN0 = 10^(EbN0_dB(i)/10);  % Conversão de dB para linear
    EsN0 = EbN0 * k;  % Razão de energia do símbolo para densidade espectral de potência de ruído
    
    sigma = sqrt(45/(2*EsN0));  % Desvio padrão do ruído
    disp(sigma)
    seq = randsample(simbolos, N, true);  % Símbolos aleatórios
    
    %s = 2*seq - (M-1);  % Modulação 4-PAM
    ruido = sigma * randn(1, N);  % Ruído AWGN
    r = seq + ruido;  % Sinal recebido
    
    demod = simbolos(arrayfun(@(x) find(abs(simbolos-x) == min(abs(simbolos-x)), 1),r));  % Demodulação
    %symbols_hat = min(max(symbols_hat, 0), M-1);  % Limitar ao intervalo de símbolos válidos
    N_erro=sum(seq ~= demod);
    Pe_sim(i) = N_erro / N;  % SER
end

% SER teórico
Pe_theory = 2*(M-1)/M.*qfunc(sqrt(6*k/(M^2-1)*10.^(EbN0_dB/10)));

% Plotar resultados
figure;
semilogy(EbN0_dB, Pe_sim, 'o-', 'DisplayName', 'Simulação');
hold on;
semilogy(EbN0_dB, Pe_theory, 's-', 'DisplayName', 'Teórico');
xlabel('E_b/N_0 (dB)');
ylabel('Taxa de Erro de Símbolo (SER)');
legend show;
title('Taxa de Erro de Símbolo 4-PAM');
grid on;