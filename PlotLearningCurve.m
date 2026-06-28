function PlotLearningCurve(squared_errors, excess_errors, iterations)
MSE = mean(squared_errors);
EMSE = mean(excess_errors);
figure;
plot(1:iterations, 10*log10(MSE), 'LineWidth', 1.5);
title('Curva de Aprendizado - MSE');
xlabel('Iterações Temporais (n)');
ylabel('MSE (dB)');
grid on;

figure;
plot(1:iterations, 10*log10(EMSE), 'r', 'LineWidth', 1.5);
title('Curva do Erro Médio Quadrático em Excesso (EMSE)');
xlabel('Iterações Temporais (n)');
ylabel('EMSE (dB)');
grid on;
end