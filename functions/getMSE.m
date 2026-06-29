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