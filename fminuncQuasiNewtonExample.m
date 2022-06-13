%% 
% Sofiya Makarenka
% Ada Kawała

%ogólny koncept naszego rozwiązania znaliezłyśmy na forumie matlabowym: 
%https://www.mathworks.com/matlabcentral/answers/658023-how-to-plot-the-trajectory-that-fminsearch-follows

%METODA QUASI_NEWTON

%%
clc
clear 
close all

%%
numberOfStartPoints = 4;
theoreticalResult = [0,0];
error = zeros(1,4);
realResults = zeros(2,4);
iterationsVector = zeros(1,4);

numberOfPlots = 1;
for k = 1:numberOfStartPoints
    [trajectory, fval, iterations, realResult]=doIt(k);
    %%
    %zapisanie końcowego wyniku w wektorze wyników
    realResults(1,k) = realResult(1,1);
    realResults(2,k) = realResult(1,2);
    %% 
    %liczba iteracji dla każdego puntku startowego
    iterationsVector(1,k) = iterations;
    %%
    %Dla każdego punktu startowego, odrębnie, wykres 2D konturowy funkcji z
    %naniesionymi trajektoriami 
    figure(numberOfPlots)
    fcontour(@rosenbrock, [-2.5 3 -3 3],'MeshDensity',50,...
                'LineWidth', 2, 'LevelList', 1:25:300);     
    hold on
    
    for i = 1 : size(trajectory, 2) - 1
        x = trajectory(1, i);
        y = trajectory(2, i);
        u = trajectory(1, i + 1) - trajectory(1, i);
        v = trajectory(2, i + 1) - trajectory(2, i);
        quiver(x, y, u, v, 0,'black', 'filled', 'linewidth', 1, 'MaxHeadSize', 2);
        scatter(trajectory(1, i), trajectory(2, i), 5, 'black', 'filled');
    end
    
    %kropka zielona - punkt startowy
    scatter(trajectory(1, 1), trajectory(2, 1), 35, 'green', 'filled');
    %kropka czerwona - rozwiązanie znalezione przez algorytm
    scatter(trajectory(1, end), trajectory(2, end), 35, 'red', 'filled');

    hold off
    xlabel('x');
    ylabel('y');
    title 'Rosenbrock solution via fminunc quasi-newton with Estimated Derivatives'
    
    numberOfPlots = numberOfPlots + 1;
    %%
    %logarytmeczna zależność f(x,y) od numeru iteracji
    figure(numberOfPlots)
    semilogy(1:iterations+1, fval, '*-');
    title 'Function value at each iteration'
    xlabel('Iteration');
    ylabel('Function Value');   
    
    %%
    %obliczenie błędu średniokwadratowego
    error(1,k) = (realResult(1,1) - theoreticalResult(1,1))^2 +(realResult(1,2) - theoreticalResult(1,2))^2;
    numberOfPlots = numberOfPlots + 1;
end 

%%
%funkcja dzięki, której zbieramy informację o poszczególnych krokach
%algorytmu
function [historyXY, historyF, iterations, result]=doIt(k)

options = optimoptions('fminunc','Display','iter',...
    'outputFcn',@out,'Algorithm','quasi-newton');
%punkty startowe zostałe umieszczone w jednym wektorze
X = [1, 0, -2, -2]; %X1 = 1, X2 = 0 i t.d..
Y = [0, -2, -2, 0]; %Y1 = 0, Y2 = -2 i t.d..

historyXY=[];
historyF=[];

[x,fval,eflag,output] = fminunc(@rosenbrock_wrapper, [X(k), Y(k)], options);
iterations = output.iterations;
result = x;
 
    %Output Function
    %nasza funkcja outputowa, służąca po to, żeby mieć dostęp do poszczególnych
    %kroków iteracji(wartości [x,y])
    %rozwiązanie powstało dzięki pliku matlabowemu - bananoout.m
    function stop = out(x, ~, state)
        stop = false;
        switch state
            case 'init'
            case 'interrupt'
            case 'iter'
                historyXY = [historyXY,x(:)];
                historyF =[ historyF, rosenbrock(x(1), x(2))];
            case 'done'
        
        end
    end
end

%%
%Rosenbrock Function - funkcja bananowa
function f = rosenbrock(x, y)
    a = -1; 
    b = -1;
    f = (1 - x  + a) ^ 2 + 100 * (y - b - (x - a) ^ 2) ^ 2;
end

%%
%Rosenbrock Wrapper
function f = rosenbrock_wrapper(X)
    f = rosenbrock(X(:, 1), X(:, 2));
end
