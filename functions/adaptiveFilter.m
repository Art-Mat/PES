function k = adaptiveFilter(x, y, m)
addpath("constants/")

% alocação de matriz para salvar os erros 
squared_errors = zeros(PartOneConsts.REALIZATIONS, PartOneConsts.ITERATIONS);

for i = 1:size(x, 1) % Loop das realizações
    u = zeros(m, 1); % Entrada do sistema
    v = zeros(m, 1); % Pesos preditos
    for j = 1:size(x, 2)
        % Atualizando o vetor de entrada
        u(2:end) = u(1:end-1);
        u(1) = x(i, j);
        d = y(i, j);
        
        error = d - dot(u, v);
        squared_errors(i, j) = error^2; % Erro quadrático da predição considerando ruído da saída
        v = v + PartTwoConsts.MU * u * error; % Atualiza os pesos preditos
    end
end
k = v; % Retorna os pesos preditos após todas as iterações
PlotMSECurve(squared_errors, size(x, 2))
end