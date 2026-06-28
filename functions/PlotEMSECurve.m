function PlotEMSECurve(excess_errors, iterations)
EMSE = mean(excess_errors);
figure;
plot(1:iterations, 10*log10(EMSE), 'r', 'LineWidth', 1.5);
title('Curva do Erro Médio Quadrático em Excesso (EMSE)');
xlabel('Iterações Temporais (n)');
ylabel('EMSE (dB)');
grid on;
end