close all; 
clear all;

load ber_COFDM_SOVA.log;
% load SNR_to_RBIR.txt;

figure(1);
set(gca, 'fontsize', 14);
semilogy(0, 1, 'w');
hold on;
grid on;
semilogy(0, 1, 'w');
semilogy(0, 1, 'w');
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,2), '-bo',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,3), '-go',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,4), '-ko',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,5), '-go',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,6), '-r^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,7), '-k^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,8), '-rs',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,9), '-m^',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,10), '-ms',  'LineWidth', 2.0, 'MarkerSIze', 10);
semilogy(ber_COFDM_SOVA(:,1), ber_COFDM_SOVA(:,11), '-bs',  'LineWidth', 2.0, 'MarkerSIze', 10);

% r=0:1:10;
% Pb=0.5.*erfc(sqrt(10.^(r./10)));
% semilogy(r, Pb, '-.kx', 'LineWidth', 2.0, 'MarkerSIze', 10);

% load ber_non_tail_biting.log;
% semilogy(ber_non_tail_biting(:,1), ber_non_tail_biting(:,2), '-.k', 'LineWidth', 2.0, 'MarkerSIze', 10);

title('BLER Performance of Coded-OFDM in ITU VA Channel Snaps');
xlabel('E_b/N_0 (dB)');
ylabel('Block Error Rate');
% axis([0 10 1e-6 1]);

legend('\{171,133\}', 'SOVA Decoder', 'FFT Size = 1024', 3);
%print -djpeg100 COFDM_SOVA_ituva_V60.jpg;
