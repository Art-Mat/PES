function k = knownPlant()
rng(PartOneConsts.SEED);
w_0 = [-0.5247, -0.2060, 0.3324, -0.2631, 0.2358, 0.0304, 0.3525, -0.4812, -0.0485, -0.2964]; % Pesos da planta
ordem = length(w_0);
squared_errors = zeros(PartOneConsts.REALIZATIONS, PartOneConsts.ITERATIONS);
for i = 1:PartOneConsts.REALIZATIONS % Loop das realizações
    u = zeros(ordem, 1); % Entrada do sistema
    v = zeros(ordem, 1); % Pesos preditos
    for j = 1:PartOneConsts.ITERATIONS
        % Atualizando o vetor de entrada
        u = circshift(u, 1);
        u(1) = randn();

        d = dot(u, w_0); % Vetor de saída da planta
        error = d - dot(u, v); % Erro de predição
        v = v + PartOneConsts.MU * u * error; % Atualiza os pesos preditos
        squared_errors(i, j) = error^2;
    end
end
k = v; % Retorna os pesos preditos após todas as iterações
PlotLearningCurve(squared_errors, PartOneConsts.ITERATIONS)
end