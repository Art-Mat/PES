function k = knownPlant()
addpath("constants/")
rng(PartOneConsts.SEED);
w_0 = [-0.5247, -0.2060, 0.3324, -0.2631, 0.2358, 0.0304, 0.3525, -0.4812, -0.0485, -0.2964]; % Pesos da planta
ordem = length(w_0); % ordem da planta

% alocação de matriz para salvar os erros 
squared_errors = zeros(PartOneConsts.REALIZATIONS, PartOneConsts.ITERATIONS);
excess_errors = zeros(PartOneConsts.REALIZATIONS, PartOneConsts.ITERATIONS);

for i = 1:PartOneConsts.REALIZATIONS % Loop das realizações
    u = zeros(ordem, 1); % Entrada do sistema
    v = zeros(ordem, 1); % Pesos preditos
    for j = 1:PartOneConsts.ITERATIONS
        % Atualizando o vetor de entrada
        u = circshift(u, 1);
        u(1) = randn(); % Considerando um ruído de variância unitária, podemos desconsiderar no cálculo

        d = dot(u, w_0) + sqrt(PartOneConsts.VAR) * randn(); % Vetor de saída da planta
        error = d - dot(u, v);
        excess_errors(i, j) = (dot(u, w_0) - dot(u, v))^2; % Erro da predição desconsiderando ruído de saída
        squared_errors(i, j) = error^2; % Erro quadrático da predição considerando ruído da saída
        v = v + PartOneConsts.MU * u * error; % Atualiza os pesos preditos
    end
end
k = v; % Retorna os pesos preditos após todas as iterações
PlotMSECurve(squared_errors, PartOneConsts.ITERATIONS)
PlotEMSECurve(excess_errors, PartOneConsts.ITERATIONS)
end