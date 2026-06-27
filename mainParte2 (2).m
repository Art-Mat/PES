function mainProjeto()
clear; close all; clc;

scriptDir = fileparts(mfilename('fullpath'));
if isempty(scriptDir)
    scriptDir = pwd;
end

outputDir = fullfile(scriptDir, 'figuras');
if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

fprintf('PARTE 1\n');
runParte1(scriptDir, outputDir);

fprintf('PARTE 2\n');
runParte2(scriptDir, outputDir);

end

% PARTE 1
function runParte1(scriptDir, outputDir)

% 1. Planta e parametros
w0 = [-0.5247; -0.2060; 0.3324; -0.2631; 0.2358; ...
       0.0304;  0.3525; -0.4812; -0.0485; -0.2964];

L = numel(w0);
M = 1000;
N = 1500;
mu = 0.01;
sigma_v2 = 1e-3;
rng(3431, 'twister');

% 2. Leitura do arquivo com entrada branca
files = dir(fullfile(scriptDir, 'whites*.mat'));
assert(~isempty(files), ...
    'Nenhum arquivo whites*.mat foi encontrado na pasta atual.');

S = load(fullfile(files(1).folder, files(1).name), 'X');
assert(isfield(S, 'X'), 'O arquivo deve conter a variavel X.');

X = double(S.X);
[nRows, nCols] = size(X);

assert(nRows * nCols >= M * N, ...
    'O arquivo X nao contem amostras suficientes.');

nBlocksPerRow = floor(nCols / N);
u = zeros(nRows * nBlocksPerRow, N);

idxReal = 1;
for r = 1:nRows
    for b = 1:nBlocksPerRow
        cols = (b - 1) * N + (1:N);
        u(idxReal, :) = X(r, cols);
        idxReal = idxReal + 1;
    end
end

u = u(1:M, :);

mean_u = mean(u, 'all');
sigma_u2 = mean(u(:).^2);

maxLag = 10;
r_u = zeros(maxLag + 1, 1);
for lag = 0:maxLag
    r_u(lag + 1) = mean(u(:, 1:N-lag) .* u(:, 1+lag:N), 'all');
end

fprintf('Media da entrada: %.6e\n', mean_u);
fprintf('Potencia/variancia da entrada: %.8f\n', sigma_u2);
fprintf('Autocorrelacoes estimadas r_u[0:10]:\n');
disp(r_u.');

% 3. Modelo da planta e inicializacao do LMS
v = sqrt(sigma_v2) * randn(M, N);

W = zeros(L, M);
reg = zeros(L, M);

mse = zeros(1, N);
emse = zeros(1, N);
msd = zeros(1, N);
wMean = zeros(L, N);

% 4. Recursao LMS
for n = 1:N
    reg(2:end, :) = reg(1:end-1, :);
    reg(1, :) = u(:, n).';

    d = w0.' * reg + v(:, n).';
    y = sum(W .* reg, 1);
    e = d - y;

    Werr = repmat(w0, 1, M) - W;
    ea = sum(Werr .* reg, 1);

    mse(n) = mean(e.^2);
    emse(n) = mean(ea.^2);
    msd(n) = mean(sum(Werr.^2, 1));
    wMean(:, n) = mean(W, 2);

    W = W + mu * reg .* repmat(e, L, 1);
end

% 5. Resultados numericos
wFinal = mean(W, 2);
wStd = std(W, 0, 2);

muMaxMS = 2 / ((L + 2) * sigma_u2);

emseInfTheory = ...
    mu * L * sigma_u2 * sigma_v2 / ...
    (2 - mu * (L + 2) * sigma_u2);

mseInfTheory = sigma_v2 + emseInfTheory;
msdInfTheory = emseInfTheory / sigma_u2;

steadyIdx = (N - 299):N;
mseInfSim = mean(mse(steadyIdx));
emseInfSim = mean(emse(steadyIdx));
msdInfSim = mean(msd(steadyIdx));

fprintf('\nLimite conservador de estabilidade em media quadratica: mu < %.6f\n', ...
    muMaxMS);
fprintf('MSE de regime (simulado):  %.8e\n', mseInfSim);
fprintf('MSE de regime (teorico):   %.8e\n', mseInfTheory);
fprintf('EMSE de regime (simulado): %.8e\n', emseInfSim);
fprintf('EMSE de regime (teorico):  %.8e\n', emseInfTheory);
fprintf('MSD de regime (simulado):  %.8e\n\n', msdInfSim);

T = table((0:L-1).', w0, wFinal, wFinal-w0, wStd, ...
    'VariableNames', {'k','wVerdadeiro','wMedioFinal', ...
                      'ErroMedioFinal','DesvioPadrao'});
disp(T);

writetable(T, fullfile(outputDir, 'parte1ResultadosCoeficientes_matlab.csv'));

% 6. Graficos
smoothWindow = 10;
msePlot = movmean(mse, smoothWindow);
emsePlot = movmean(emse, smoothWindow);
msdPlot = movmean(msd, smoothWindow);

figure;
plot(0:N-1, 10*log10(max(msePlot, eps)), 'LineWidth', 1.2);
hold on;
plot(0:N-1, 10*log10(max(emsePlot, eps)), 'LineWidth', 1.2);
yline(10*log10(mseInfTheory), '--', 'MSE teorico');
yline(10*log10(emseInfTheory), '--', 'EMSE teorico');
grid on;
xlabel('Iteracao n');
ylabel('Potencia do erro (dB)');
title('MSE e EMSE - identificação LMS');
legend('MSE', 'EMSE', 'Location', 'best');
exportgraphics(gcf, fullfile(outputDir, 'parte1MseEmseMatlab.png'), 'Resolution', 180);

figure;
plot(0:N-1, wMean, 'LineWidth', 1.0);
hold on;
for k = 1:L
    yline(w0(k), '--', 'HandleVisibility', 'off');
end
grid on;
xlabel('Iteracao n');
ylabel('Coeficiente');
title('Convergência dos coeficientes médios');
legend(compose('w_%d', 0:L-1), 'Location', 'eastoutside');
exportgraphics(gcf, fullfile(outputDir, 'parte1CoeficientesMatlab.png'), 'Resolution', 180);

figure;
plot(0:N-1, 10*log10(max(msdPlot, eps)), 'LineWidth', 1.2);
hold on;
yline(10*log10(msdInfTheory), '--', 'MSD teorico');
grid on;
xlabel('Iteracao n');
ylabel('MSD (dB)');
title('Desvio quadrático médio dos coeficientes');
exportgraphics(gcf, fullfile(outputDir, 'parte1MsdMatlab.png'), 'Resolution', 180);

end

% PARTE 2
function runParte2(scriptDir, outputDir)

cfg.dataDir = scriptDir;
cfg.outputDir = outputDir;

cfg.layout = 'auto';
cfg.maxCoeffs = 80;
cfg.orderTolerance = 0.02;
cfg.ridge = 1e-8;

cfg.muWhite = 0.35;
cfg.muCorrel = 0.08;
cfg.nlmsDelta = 1e-6;
cfg.learningMaxSamples = inf;
cfg.smoothWindow = 200;

cases = {
    'whites.mat', 'branco', cfg.muWhite;
    'correl.mat', 'correlacionado', cfg.muCorrel
};

results = struct([]);

for k = 1:size(cases, 1)
    matFile = fullfile(cfg.dataDir, cases{k, 1});
    label = cases{k, 2};
    mu = cases{k, 3};

    results(k).label = label;
    results(k).file = matFile;
    results(k).mu = mu;
    results(k).case = analyzeCase(matFile, label, mu, cfg);
end

plotComparison(results, cfg);

fprintf('\nResumo Parte 2:\n');
for k = 1:numel(results)
    r = results(k).case;
    fprintf('Caso %-15s | layout %-11s | L = %d coef. | ordem = %d | NMSE val final = %.4g\n', ...
        results(k).label, r.layoutName, r.L, r.order, r.valNMSE(r.L));
    fprintf('w_hat = [');
    fprintf(' %.8g', r.wHat);
    fprintf(' ]^T\n\n');
end

end

function out = analyzeCase(matFile, label, mu, cfg)
[U, D, layoutName] = loadSignals(matFile, cfg);

numRuns = size(U, 1);
numSamples = min(size(U, 2), size(D, 2));
U = U(:, 1:numSamples);
D = D(:, 1:numSamples);

trainRows = 1:max(1, floor(0.70 * numRuns));
valRows = (max(trainRows) + 1):numRuns;
if isempty(valRows)
    valRows = trainRows;
end

[L, order, trainNMSE, valNMSE, allW] = chooseOrder(U, D, trainRows, valRows, cfg);
wHat = allW{L};

fprintf('CASO %s: Numero de coeficientes escolhido L = %d\n', upper(label), L);
fprintf('CASO %s: Ordem estimada da planta = L - 1 = %d\n', upper(label), order);

[mseLearning, wNlmsMean] = nlmsLearningCurve(U, D, L, mu, cfg);

out.U = U;
out.D = D;
out.layoutName = layoutName;
out.L = L;
out.order = order;
out.wHat = wHat;
out.wNlmsMean = wNlmsMean;
out.trainNMSE = trainNMSE;
out.valNMSE = valNMSE;
out.mseLearning = mseLearning;

plotOrderCurve(trainNMSE, valNMSE, L, label, cfg);
plotLearningCurve(mseLearning, label, cfg);
plotPlantVector(wHat, label, cfg);

end

function [U, D, layoutName] = loadSignals(matFile, cfg)
if ~exist(matFile, 'file')
    error('Arquivo nao encontrado: %s', matFile);
end

S = load(matFile);
names = fieldnames(S);

if isfield(S, 'X') && isfield(S, 'Y')
    U = double(S.X);
    D = double(S.Y);
    layoutName = 'X=entrada, Y=saida';
    return;
end

if isfield(S, 'X')
    X = S.X;
else
    X = S.(names{1});
end

X = double(X);
[U, D, layoutName] = splitInputOutput(X, cfg);

end

function [U, D, chosenName] = splitInputOutput(X, cfg)
candidates = {};

if mod(size(X, 2), 2) == 0
    n = size(X, 2);
    candidates{end + 1} = makeCandidate('halves', X(:, 1:n/2), X(:, n/2+1:end));
    candidates{end + 1} = makeCandidate('interleaved', X(:, 1:2:end), X(:, 2:2:end));
end

if mod(size(X, 1), 2) == 0
    candidates{end + 1} = makeCandidate('rowpairs', X(1:2:end, :), X(2:2:end, :));
end

if isempty(candidates)
    error('Nao foi possivel montar candidatos de entrada/saida a partir de X.');
end

if ~strcmpi(cfg.layout, 'auto')
    idx = find(strcmpi(cfg.layout, cellfun(@(c) c.name, candidates, 'UniformOutput', false)), 1);
    if isempty(idx)
        error('Layout cfg.layout = %s nao existe para esta matriz.', cfg.layout);
    end
else
    nmse = zeros(numel(candidates), 1);
    for k = 1:numel(candidates)
        maxL = min(20, cfg.maxCoeffs);
        rows = 1:min(20, size(candidates{k}.U, 1));
        [~, vNMSE] = quickFitScore(candidates{k}.U, candidates{k}.D, rows, maxL, cfg);
        nmse(k) = min(vNMSE);
    end
    [~, idx] = min(nmse);
end

U = candidates{idx}.U;
D = candidates{idx}.D;
chosenName = candidates{idx}.name;

end

function c = makeCandidate(name, U, D)
n = min(size(U, 2), size(D, 2));
c.name = name;
c.U = double(U(:, 1:n));
c.D = double(D(:, 1:n));
end

function [bestL, valNMSE] = quickFitScore(U, D, rows, maxL, cfg)
numRuns = numel(rows);
trainRows = rows(1:max(1, floor(0.7 * numRuns)));
valRows = rows(max(trainRows)+1:end);

if isempty(valRows)
    valRows = trainRows;
end

valNMSE = inf(maxL, 1);

for L = 1:maxL
    w = estimateWiener(U, D, trainRows, L, cfg);
    [~, valNMSE(L)] = evaluateModel(U, D, valRows, w);
end

[~, bestL] = min(valNMSE);

end

function [LChosen, orderChosen, trainNMSE, valNMSE, allW] = chooseOrder(U, D, trainRows, valRows, cfg)
maxL = min(cfg.maxCoeffs, min(size(U, 2), size(D, 2)) - 1);
trainNMSE = inf(maxL, 1);
valNMSE = inf(maxL, 1);
allW = cell(maxL, 1);

for L = 1:maxL
    w = estimateWiener(U, D, trainRows, L, cfg);
    allW{L} = w;
    [~, trainNMSE(L)] = evaluateModel(U, D, trainRows, w);
    [~, valNMSE(L)] = evaluateModel(U, D, valRows, w);
end

minNMSE = min(valNMSE);
threshold = minNMSE * (1 + cfg.orderTolerance);
LChosen = find(valNMSE <= threshold, 1, 'first');
orderChosen = LChosen - 1;

end

function w = estimateWiener(U, D, rows, L, cfg)
[ruu, rdu] = estimateCorrelations(U, D, rows, L);
R = toeplitz(ruu(1:L));
p = rdu(1:L);

lambda = cfg.ridge * max(trace(R) / L, eps);
w = (R + lambda * eye(L)) \ p(:);

end

function [ruu, rdu] = estimateCorrelations(U, D, rows, L)
numSamples = min(size(U, 2), size(D, 2));
startIdx = L;

ruu = zeros(L, 1);
rdu = zeros(L, 1);
count = 0;

for rr = rows
    u = U(rr, :);
    d = D(rr, :);
    idx = startIdx:numSamples;
    count = count + numel(idx);

    for lag = 0:(L - 1)
        ruu(lag + 1) = ruu(lag + 1) + sum(u(idx) .* u(idx - lag));
        rdu(lag + 1) = rdu(lag + 1) + sum(d(idx) .* u(idx - lag));
    end
end

ruu = ruu / max(count, 1);
rdu = rdu / max(count, 1);

end

function [mse, nmse] = evaluateModel(U, D, rows, w)
L = numel(w);
numSamples = min(size(U, 2), size(D, 2));
sse = 0;
sd = 0;
count = 0;

for rr = rows
    u = U(rr, 1:numSamples);
    d = D(rr, 1:numSamples);
    yHat = filter(w(:).', 1, u);
    e = d(L:end) - yHat(L:end);
    sse = sse + sum(e.^2);
    sd = sd + sum(d(L:end).^2);
    count = count + numel(e);
end

mse = sse / max(count, 1);
nmse = sse / max(sd, eps);

end

function [mseLearning, wMean] = nlmsLearningCurve(U, D, L, mu, cfg)
numRuns = size(U, 1);
numSamples = min(size(U, 2), size(D, 2));
numSamples = min(numSamples, cfg.learningMaxSamples);

mseLearning = zeros(numSamples, 1);
wSum = zeros(L, 1);

for rr = 1:numRuns
    w = zeros(L, 1);
    u = U(rr, 1:numSamples);
    d = D(rr, 1:numSamples);

    for n = L:numSamples
        x = u(n:-1:n-L+1).';
        e = d(n) - w.' * x;
        w = w + (mu / (cfg.nlmsDelta + x.' * x)) * x * e;
        mseLearning(n) = mseLearning(n) + e.^2;
    end

    wSum = wSum + w;
end

mseLearning = mseLearning / numRuns;
mseLearning(1:L-1) = nan;
wMean = wSum / numRuns;

end

function plotOrderCurve(trainNMSE, valNMSE, LChosen, label, cfg)
fig = figure('Color', 'w');
Lvals = 1:numel(valNMSE);

semilogy(Lvals, trainNMSE, 'LineWidth', 1.4);
hold on;
semilogy(Lvals, valNMSE, 'LineWidth', 1.4);

yl = ylim;
plot([LChosen LChosen], yl, '--k', 'LineWidth', 1.2);
ylim(yl);

grid on;
xlabel('Numero de coeficientes L');
ylabel('NMSE');
title(sprintf('Escolha da ordem - caso %s', label), 'Interpreter', 'none');
legend('Treino', 'Validacao', 'L escolhido', 'Location', 'best');

saveFigure(fig, cfg.outputDir, sprintf('parte2ordem%s', label));

end

function plotLearningCurve(mseLearning, label, cfg)
fig = figure('Color', 'w');
valid = ~isnan(mseLearning);
mseSmooth = mseLearning;

if cfg.smoothWindow > 1
    mseSmooth(valid) = movmean(mseLearning(valid), cfg.smoothWindow);
end

semilogy(mseLearning, 'Color', [0.75 0.75 0.75], 'LineWidth', 0.8);
hold on;
semilogy(mseSmooth, 'b', 'LineWidth', 1.5);

grid on;
xlabel('Iteracao i');
ylabel('MSE(i)');
title(sprintf('Curva de aprendizado - caso %s', label), 'Interpreter', 'none');
legend('MSE instantaneo medio', 'MSE suavizado', 'Location', 'best');

saveFigure(fig, cfg.outputDir, sprintf('parte2MseAprendizado%s', label));

end

function plotPlantVector(wHat, label, cfg)
fig = figure('Color', 'w');

stem(0:numel(wHat)-1, wHat, 'filled', 'LineWidth', 1.2);

grid on;
xlabel('Indice do coeficiente');
ylabel('w estimado');
title(sprintf('Vetor estimado da planta - caso %s', label), 'Interpreter', 'none');

saveFigure(fig, cfg.outputDir, sprintf('parte2VetorEstimado%s', label));

end

function plotComparison(results, cfg)
fig = figure('Color', 'w');

for k = 1:numel(results)
    mse = results(k).case.mseLearning;
    valid = ~isnan(mse);
    mseSmooth = mse;
    mseSmooth(valid) = movmean(mse(valid), cfg.smoothWindow);
    semilogy(mseSmooth, 'LineWidth', 1.4);
    hold on;
end

grid on;
xlabel('Iteracao i');
ylabel('MSE(i) suavizado');
title('Comparação das curvas de aprendizado');
legend({results.label}, 'Location', 'best');

saveFigure(fig, cfg.outputDir, 'parte2ComparacaoMSE');

end

function saveFigure(fig, outputDir, baseName)
pdfFile = fullfile(outputDir, [baseName '.pdf']);
pngFile = fullfile(outputDir, [baseName '.png']);

try
    exportgraphics(fig, pdfFile, 'ContentType', 'vector');
    exportgraphics(fig, pngFile, 'Resolution', 300);
catch
    set(fig, 'PaperPositionMode', 'auto');
    print(fig, pdfFile, '-dpdf', '-bestfit');
    print(fig, pngFile, '-dpng', '-r300');
end

end