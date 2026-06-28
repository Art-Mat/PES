classdef PartOneConsts
    properties (Constant = true)
        REALIZATIONS = 1000; % Número de realizações (definido no enunciado)
        ITERATIONS = 500; % Número de iterações até convergência
        SEED = 0; % Seed do motor de número aleatório
        MU = 0.02; % Passo de adaptação do algoritmo
        VAR = 0.001 % Variância do ruído da saída
    end
end