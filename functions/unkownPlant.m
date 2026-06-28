function k = unkownPlant()
addpath("data")
addpath("functions/")
x = load("whites.mat").X;
y = load("whites.mat").Y;

getOrder(x, y) % Gráfico para questão a)


end

function getOrder(x, y)
% Plota o MSE da Solução de Wiener para diferentes ordens
    MSE = zeros(PartTwoConsts.MAX_ORDER, 1);
    for i = 1:PartTwoConsts.MAX_ORDER
            % Inicialização das matrizes e vetores usadas na iteração
            corr_matrix = zeros(i, i);
            corr_vec = zeros(i ,1);
            for j = 1:size(x, 1)
                u = zeros(i, 1);
                for k = 1:size(x, 2)
                    % Adicionando o novo valor aos vetores entrada e saída
                    u(2:end) = u(1:end-1);
                    u(1) = x(j, k);
                    d = y(j, k);
        
                    % Realizando o acúmulo para a correlação entre entradas e
                    % entrada/saída
                    corr_matrix = corr_matrix + u * u';
                    corr_vec = corr_vec + u * d;
                end
            end
            corr_matrix = corr_matrix / (size(x, 1) * size(x, 2));
            corr_vec = corr_vec / (size(x, 1) * size(x, 2));
    
            w = corr_matrix \ corr_vec; % Solução de wiener para o sistema
            MSE(i) = getMSE(w, x, y);
    end
    figure;
    scatter(0:PartTwoConsts.MAX_ORDER - 1, 10*log10(MSE), 'LineWidth', 1.5);
    title('MSE x Ordem');
    xlabel('Ordem do sistema (M)');
    ylabel('MSE (dB)');
    grid on;
    
end

function MSE = getMSE(w, X, Y)
    sce = 0; % Square cumulative error
    
    for r = 1:size(X, 1) % Loop das 100 realizações
        u = zeros(length(w), 1);
        for k = 1:size(X, 2) % Loop temporal
            u(2:end) = u(1:end-1);
            u(1) = X(r, k);
            % Acumula o erro quadrático
            sce = sce + (Y(r, k) - w' * u)^2;
        end
    end
    
    % Divide pelo número total de amostras processadas
    MSE = sce / (size(X, 1) * size(X, 2));
end