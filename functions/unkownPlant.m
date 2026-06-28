function k = unkownPlant()
addpath("data")
addpath("functions/")
load("whites.mat")
% alocação de matriz para salvar os erros 
squared_errors = zeros(PartTwoConsts.REALIZATIONS, PartTwoConsts.ITERATIONS);

k = v; % Retorna os pesos preditos após todas as iterações
PlotLearningCurve(squared_errors, excess_errors, PartOneConsts.ITERATIONS)
end

function M = getOrder(x, y)
    for i = 1:PartTwoConsts.MAX_ORDER
        corr_matrix = zeros(m, n);
        corr_vec = zeros(m ,1);
        u = zeros(i);
        d = zeros(i);
        for j = 1:length(x);
            % Adicionando o novo valor aos vetores entrada e saída
            u = circshift(u, 1);
            u(1) = x(1, j);
            d = circshift(d, 1);
            d(1) = y(1, j);

            % Realizando o acúmolo para a correlação entre entradas e
            % entrada/saída
            corr_matrix = corr_matrix + u' * u;
            corr_vec = corr_vec + u' * d;
        end
        corr_matrix = corr_matrix / length(x);
        corr_vector = corr_vector / length(x);

        w = 
    end
end