%% Clean up
clc; clear variables;

%% Mathematical Functions
f = @(x) 1/(1+(1/x)^4);
lx = @(x) x*log(x);
lx_sol = log(4)-3/4;

%% Naming scheme: <value type>_<function type><n where 10^n iterations>
% 
% value types:
% I  -> Integral
% E  -> Error
% S  -> Slope
%
% function types:
% ls -> Left riemann Sum
% rs -> Right riemann Sum
% ms -> Midpoint riemann Sum
% ts -> Trapezoid method (Sum)
% ss -> Simpsons rule (Sum)
%%% %%% %%% %%% %%% %%% %%%

%% Integral Approximation

size = 6; % You can manipulate the variable to extend the graph
          % Warning: It takes noticibly longer for the program 
          %          to run on higher 'size's.
functionNames = [
    "left Riemann sum";
    "right Riemann sum";
    "middle Riemann";
    "trapezoid method";
    "Simpsons rule"
    ];
numOfFunctions = numel(functionNames);
operation = [
    "approximation";
    "error";
    ];

I_ls = zeros(1, size);
I_rs = zeros(1, size);
I_ms = zeros(1, size);
I_ts = zeros(1, size);
I_ss = zeros(1, size);

for i = 1:size
    I_ls(1, i) = leftRiemann(lx, 1, 2, 10^i);
    I_rs(1, i) = rightRiemann(lx, 1, 2, 10^i);
    I_ms(1, i) = middleRiemann(lx, 1, 2, 10^i);
    I_ts(1, i) = trapezoidMethod(lx, 1, 2, 10^i);
    I_ss(1, i) = simpsonsRule(lx, 1, 2, 10^i);
end

I    = [I_ls;
        I_rs;
        I_ms;
        I_ts;
        I_ss];

printData(functionNames, numOfFunctions, I, size, operation(1));


%% Error Calculation
E_ls = zeros(1, size);
E_rs = zeros(1, size);
E_ms = zeros(1, size);
E_ts = zeros(1, size);
E_ss = zeros(1, size);
X    = zeros(1, size);

for i = 1:size
    X(i)    = 10^i;
    E_ls(i) = abs(lx_sol - I_ls(i))/lx_sol;
    E_rs(i) = abs(lx_sol - I_rs(i))/lx_sol;
    E_ms(i) = abs(lx_sol - I_ms(i))/lx_sol;
    E_ts(i) = abs(lx_sol - I_ts(i))/lx_sol;
    E_ss(i) = abs(lx_sol - I_ss(i))/lx_sol;
end

E    = [E_ls;
        E_rs;
        E_ms;
        E_ts;
        E_ss];

printData(functionNames, numOfFunctions, E, size, operation(2));


%% Plotting Error
figure; %grid on; box on; hold on;

loglog(E_ls, X,  ...
       E_rs, X,  ...
       E_ms, X,  ...
       E_ts, X,  ...
       E_ss, X   );
hold on;

legend( "Left Riemann Error", "Right Riemann Error", ...
    "Middle Riemann Error", "Trapezoid Method", "Simpson's Rule" );

xlabel("Error Value");
ylabel("Number of Iterations");
title("Error Values v. Number of Iterations");


%% Slope Calculation
S_ls = log10(E_ls(1)/E_ls(size)) / log10(X(size)/X(1));
S_rs = log10(E_rs(1)/E_rs(size)) / log10(X(size)/X(1));
S_ms = log10(E_ms(1)/E_ms(size)) / log10(X(size)/X(1));
S_ts = log10(E_ts(1)/E_ts(size)) / log10(X(size)/X(1));
S_ss = log10(E_ss(1)/E_ss(size)) / log10(X(size)/X(1));
S    = [S_ls; S_rs; S_ms; S_ts; S_ss];

fprintf("Slope values:\n")
for t = 1:numOfFunctions
    fprintf("The slope of the %s error is: %f\n", functionNames(t), S(t));
end


%% Print Function
function printData(nameList, nameListSize, data, dataLength, op)
    fprintf("All %ss:\n", op);
    for f = 1:nameListSize
        for s = 1:dataLength
            fprintf("The %s for 10^%d points with the %s is: %f\n", op, s, nameList(f), data(f, s));
        end
        fprintf("\n");
    end
    fprintf("\n");
end


%% Summation Functions
function leftSum = leftRiemann(f, s, e, N)
    leftSum = 0;
    h = (e-s)/N;
    for it = 0:N-1
        leftSum = leftSum + f((s+0.0)+(h*it))*h;
    end
end

function rightSum = rightRiemann(f, s, e, N)
    rightSum = 0;
    h = (e-s)/N;
    for it = 0:N-1
        rightSum = rightSum + f((s+(h/2))+(h*it))*h;
    end
end

function middleSum = middleRiemann(f, s, e, N)
    middleSum = 0;
    h = (e-s)/N;
    for it = 0:N-1
        middleSum = middleSum + f((s+h)+(h*it))*h;
    end
end

function trapezoidSum = trapezoidMethod(f, s, e, N)
    trapezoidSum = (f(s)+f(e))*(1/2);
    h = (e-s)/N;
    
    for it = 1:N-1
        trapezoidSum = trapezoidSum + f(s + (it*h));
    end

    trapezoidSum = h*trapezoidSum;
end

function simpsonSum = simpsonsRule(f, s, e, N)
    simpsonSum = 0;

    if mod(N,2) ~= 0 
        % if the gaps of Simpson's Rule are not even, this means that
        % you're going to have an outlier gap that will not be used in the
        % operation. Methods to deal with that include extrapolation or
        % deletion of data. We will go with deletion here. An extrapolation
        % would be simplified here by doing N+1.
        N = N - 1;
        fprintf("Your 'data' had an uneven number of points, so the last one was deleted to properly work with Simpson's Rule\n");
    end
    
    h = (e-s)/N; % each gap

    % otherH = h/2; % doing Simpson's rule by creating points inbetween
    % (not used)

    for it = 0:2:N-2
        % the N-2 exists here for 2 different reasons: 
        %
        % (1) a necessary shift of the data from [1, N] to [0, N-1] so the
        % calculation of the function used for Simpson's rule makes more
        % intuitive sense in that the starting point is (s+h)+(it+h) and
        % not (s+h) + ( (it-1)+h ). 
        %
        % (2) removing the last iteration of the loop prevents us from
        % operating the formula on a piece of 'data' that does not exist.
        % Despite having removed a piece of data earlier, this does not
        % change the operation here. Once you reach the end of the data,
        % having deleted the unnecessary data makes the difference between
        % evaluating 1 and 2 non existent data points

        simpsonSum = simpsonSum + (1/3)*h*f((s)+(it*h)) + (2/3)*h*f((s+h)+(it*h)) + (1/3)*h*f((s+2*h)+(it*h));

    end

end

%% Reference (not used)
%{
function generalSum = allRiemann(f, s, e, N) 
    generalSum = 0;
    h = (e-s)/N;
    for it = 0:N-1
        generalSum = generalSum + evalF(f,s,h,it);
    end
end

leftF   = @(s,h,n) f((s+0.0)+(h*n))*h; % reference (not used)
rightF  = @(s,h,n) f((s+0.5)+(h*n))*h; % reference (not used)
middleF = @(s,h,n) f((s+1.0)+(h*n))*h; % reference (not used)

%}


