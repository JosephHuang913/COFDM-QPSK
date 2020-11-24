close all; 
clear all;

load ber_COFDM_SOVA.log;
% load SNR_to_RBIR.txt;

figure(1);
set(gca, 'fontsize', 14);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,2), '-ro',  'LineWidth', 2.0, 'MarkerSIze', 10);
hold on;
grid on;
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,3), '-b^',  'LineWidth', 2.0, 'MarkerSIze', 10);
plot(0, 1, 'w');
plot(0, 1, 'w');
plot(0, 1, 'w');

% r=0:1:10;
% Pb=0.5.*erfc(sqrt(10.^(r./10)));
% semilogy(r, Pb, '-.kx', 'LineWidth', 2.0, 'MarkerSIze', 10);

% load ber_non_tail_biting.log;
% semilogy(ber_non_tail_biting(:,1), ber_non_tail_biting(:,2), '-.k', 'LineWidth', 2.0, 'MarkerSIze', 10);

title('BER and BLER Performance of Coded-OFDM in ITU Vehicular A Channel');
xlabel('E_b/N_0 (dB)');
ylabel('Error Rate');
axis([0 18 1e-5 1]);

legend('BER', 'BLER',  '\{171,133\}', 'SOVA Decoder', 'FFT Size = 1024', 3);
%print -djpeg100 COFDM_SOVA_ITUVA_V0.jpg;
