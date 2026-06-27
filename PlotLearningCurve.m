function PlotLearningCurve(squared_errors, iterations)
MSE = mean(squared_errors);
figure;
plot(1:iterations, 10*log10(MSE), 'LineWidth', 1.5);
title('Curva de Aprendizado - MSE');
xlabel('Iterações Temporais (n)');
ylabel('MSE (dB)');
grid on;
end