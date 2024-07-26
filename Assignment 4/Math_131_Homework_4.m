clc; clear variables; clear figures;

%% Anonymous Functions
f  = @(x) 10.*exp(-1.*x./5).*cos(x);
fp = @(x) -2.*exp(-1.*x./5).*(5.*sin(x)+cos(x));

s  = @(x) sin(x);
c  = @(x) cos(x);



%% Problem 1
x = 0:pi/6:pi;
y = [0, 1/2, sqrt(3)/2, 1, sqrt(3)/2, 1/2, 0];

FirstOrderFWD_List = firstOrderForward(x, y);
SecondOrderCTR_List = secondOrderCenter(x, y);
FirstOrderBWD_List = firstOrderBackward(x, y);

fprintf("\nThese are the values of x:\n");
printList(x);
fprintf("These are the values of y:\n");
printList(y);
fprintf("\n\n");

fprintf("These are the approximations:\n");
printList(FirstOrderFWD_List);
printList(SecondOrderCTR_List);
printList(FirstOrderBWD_List);
fprintf("\n\n\n");



%% Problem 2
%%% Create the 'data' and put it in the First-Order Methods
x_10  = linspace(0, 10, 10);
x_100 = linspace(0, 10, 100);
x_1000= linspace(0, 10, 1000);

y_10   = f(x_10);
y_100  = f(x_100);
y_1000 = f(x_1000);

FirstOrderFWD_List_10   = firstOrderForward(x_10, y_10);
FirstOrderFWD_List_100  = firstOrderForward(x_100, y_100);
FirstOrderFWD_List_1000 = firstOrderForward(x_1000, y_1000);

FirstOrderBWD_List_10   = firstOrderBackward(x_10, y_10);
FirstOrderBWD_List_100  = firstOrderBackward(x_100, y_100);
FirstOrderBWD_List_1000 = firstOrderBackward(x_1000, y_1000);

SecondOrderFWD_List_10   = secondOrderForward(x_10, y_10);
SecondOrderFWD_List_100  = secondOrderForward(x_100, y_100);
SecondOrderFWD_List_1000 = secondOrderForward(x_1000, y_1000);

SecondOrderCTR_List_10   = secondOrderCenter(x_10, y_10);
SecondOrderCTR_List_100  = secondOrderCenter(x_100, y_100);
SecondOrderCTR_List_1000 = secondOrderCenter(x_1000, y_1000);

SecondOrderBWD_List_10   = secondOrderBackward(x_10, y_10);
SecondOrderBWD_List_100  = secondOrderBackward(x_100, y_100);
SecondOrderBWD_List_1000 = secondOrderBackward(x_1000, y_1000);

FourthOrderCTR_List_10   = fourthOrderCenter(x_10, y_10);
FourthOrderCTR_List_100  = fourthOrderCenter(x_100, y_100);
FourthOrderCTR_List_1000 = fourthOrderCenter(x_1000, y_1000);

fprintf("Here are the new approximations:\n");
fprintf("First order forward approximations:\n");
printList(FirstOrderFWD_List_10);
printList(FirstOrderFWD_List_100);
printList(FirstOrderFWD_List_1000);
fprintf("\n");
fprintf("First order backward approximations:\n");
printList(FirstOrderBWD_List_10);
printList(FirstOrderBWD_List_100);
printList(FirstOrderBWD_List_1000);
fprintf("\n");
fprintf("Second order forward approximations:\n");
printList(SecondOrderFWD_List_10);
printList(SecondOrderFWD_List_100);
printList(SecondOrderFWD_List_1000);
fprintf("\n");
fprintf("Second order center approximations:\n");
printList(SecondOrderCTR_List_10);
printList(SecondOrderCTR_List_100);
printList(SecondOrderCTR_List_1000);
fprintf("\n");
fprintf("Second order backward approximations:\n");
printList(SecondOrderBWD_List_10);
printList(SecondOrderBWD_List_100);
printList(SecondOrderBWD_List_1000);
fprintf("\n");
fprintf("Fourth order center approximations:\n");
printList(FourthOrderCTR_List_10);
printList(FourthOrderCTR_List_100);
printList(FourthOrderCTR_List_1000);
fprintf("\n\n");
fprintf("Something to keep in mind are the true values:\n");
printList(fp(x_1000));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Error Difference Plotting %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: Sometimes the semilogy graphs do not graph in the log scale when
% plotting, and I have not been able to find what causes this inconsistency.
%
%%%%%%%%%%%%%%%%%%%%%%%%

% First Order Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure; 
subplot(2,1,1); 
semilogy(x_10(1:end-1), abs(FirstOrderFWD_List_10(1:end-1) - fp(x_10(1:end-1))), 'LineWidth', 2, 'Color', "#FF5555");
hold on;
semilogy(x_100(1:end-1), abs(FirstOrderFWD_List_100(1:end-1) - fp(x_100(1:end-1))), 'LineWidth', 2, 'Color', "#55FF55");
semilogy(x_1000(1:end-1), abs(FirstOrderFWD_List_1000(1:end-1) - fp(x_1000(1:end-1))), 'LineWidth', 2, 'Color', "#5555FF");

grid on; box on;
ylabel("f'(x) (log scale)", "FontName", "TimesNewRoman");
title("First Order Forward Approximation Errors");
FirstOrderFWDLegend = legend(...
    "First Order Approx (FWD), 10   iterations", ...
    "First Order Approx (FWD), 100  iterations", ...
    "First Order Approx (FWD), 1000 iterations", ...
    "Location", "southwest" ...
    );
FirstOrderFWDLegend.FontSize = 6;
hold off;


subplot(2,1,2);
semilogy(x_10(2:end), abs(FirstOrderBWD_List_10(2:end) - fp(x_10(2:end))), 'LineWidth', 2);
hold on;
semilogy(x_100(2:end), abs(FirstOrderBWD_List_100(2:end) - fp(x_100(2:end))), 'LineWidth', 2);
semilogy(x_1000(2:end), abs(FirstOrderBWD_List_1000(2:end) - fp(x_1000(2:end))), 'LineWidth', 2);

grid on; box on;
xlabel('x', "FontName", "TimesNewRoman");
ylabel("f '(x) (log scale)", "FontName", "TimesNewRoman");
title("First Order Backward Approximation Errors");
FirstOrderBWDLegend = legend(...
    "First Order Approx (BWD), 10   iterations", ...
    "First Order Approx (BWD), 100  iterations", ...
    "First Order Approx (BWD), 1000 iterations", ...
    "Location", "southwest" ...
    );
FirstOrderBWDLegend.FontSize = 6;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Second Order Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
subplot(2,2,1); 
semilogy(x_10(1:end-2), abs(SecondOrderFWD_List_10(1:end-2) - fp(x_10(1:end-2))), 'LineWidth', 2); 
hold on;
semilogy(x_100(1:end-2), abs(SecondOrderFWD_List_100(1:end-2) - fp(x_100(1:end-2))), 'LineWidth', 2);
semilogy(x_1000(1:end-2), abs(SecondOrderFWD_List_1000(1:end-2) - fp(x_1000(1:end-2))), 'LineWidth', 2);

grid on; box on;
ylabel("f '(x) (log scale)", "FontName", "TimesNewRoman");
title("Second Order Forward Approximation Errors", 'FontSize', 8);
SecondOrderFWDLegend = legend(...
    "Second Order Approx (FWD), 10   iterations", ...
    "Second Order Approx (FWD), 100  iterations", ...
    "Second Order Approx (FWD), 1000 iterations", ...
    "Location", "southwest" ...
    );
SecondOrderFWDLegend.FontSize = 5;
hold off;


subplot(2,2,2);
semilogy(x_10(1:end-1), abs(SecondOrderCTR_List_10(1:end-1) - fp(x_10(1:end-1))), 'LineWidth', 2);
hold on; grid on; box on;
semilogy(x_100(1:end-1), abs(SecondOrderCTR_List_100(1:end-1) - fp(x_100(1:end-1))), 'LineWidth', 2);
semilogy(x_1000(1:end-1), abs(SecondOrderCTR_List_1000(1:end-1) - fp(x_1000(1:end-1))), 'LineWidth', 2);

title("Second Order Center Approximation Errors", 'FontSize', 8);
SecondOrderCTRLegend = legend(...
    "Second Order Approx (CTR), 10   iterations", ...
    "Second Order Approx (CTR), 100  iterations", ...
    "Second Order Approx (CTR), 1000 iterations", ...
    "Location", "southwest" ...
    );
SecondOrderCTRLegend.FontSize = 5;


subplot(2,2,[3,4]);
semilogy(x_10(3:end), abs(SecondOrderBWD_List_10(3:end) - fp(x_10(3:end))), 'LineWidth', 2);
hold on; grid on; box on;
semilogy(x_100(3:end), abs(SecondOrderBWD_List_100(3:end) - fp(x_100(3:end))), 'LineWidth', 2);
semilogy(x_1000(3:end), abs(SecondOrderBWD_List_1000(3:end) - fp(x_1000(3:end))), 'LineWidth', 2);

xlabel('x', "FontName", "TimesNewRoman");
ylabel("f '(x) (log scale)", "FontName", "TimesNewRoman");
title("Second Order Backward Approximation Errors", 'FontSize', 8);
SecondOrderBWDLegend = legend(...
    "Second Order Approx (BWD), 10   iterations", ...
    "Second Order Approx (BWD), 100  iterations", ...
    "Second Order Approx (BWD), 1000 iterations", ...
    "Location", "southwest" ...
    );
SecondOrderBWDLegend.FontSize = 5;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Fourth Order Plotting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
semilogy(x_10(3:end-2), abs(FourthOrderCTR_List_10(3:end-2) - fp(x_10(3:end-2))), 'LineWidth', 2);
hold on; grid on; box on;
semilogy(x_100(3:end-2), abs(FourthOrderCTR_List_100(3:end-2) - fp(x_100(3:end-2))), 'LineWidth', 2);
semilogy(x_1000(3:end-2), abs(FourthOrderCTR_List_1000(3:end-2) - fp(x_1000(3:end-2))), 'LineWidth', 2); 

xlabel('x', "FontName", "TimesNewRoman");
ylabel("f '(x) (log scale)", "FontName", "TimesNewRoman");
title("Fourth Order Center Approximation Errors", 'FontSize', 8);
FourthOrderCTRLegend = legend(...
    "Fourth Order Approx (CTR), 10   iterations", ...
    "Fourth Order Approx (CTR), 100  iterations", ...
    "Fourth Order Approx (CTR), 1000 iterations", ...
    "Location", "southwest" ...
    );
FourthOrderCTRLegend.FontSize = 5;
hold off;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%% Print the Data
% Pretty straight forward, saves me time
function printList(list)
    for i = 1:length(list)-1
        fprintf("%10f, ", list(i));
    end
    fprintf("%10f\n", list(length(list)));
end



%% Approximation Functions
%%%%%%%%%%%%%%%%%%%%%%%%%
% I needed MATLAB to allow me to work with a flexible list, so I
% commented out the section that creates the initial list. It gives
% a warning about that, but it functions as needed.


function yp_list = firstOrderForward(x_list, y_list)
    sizeOfLists = length(x_list);  % record of the size of your lists
    %yp_list = zeros(sizeOfLists);  % make a list
    
    for n = 1:sizeOfLists-1
        h = x_list(n+1) - x_list(n);              % make h be the difference between the current x value and the next one
        yp_list(n) = (y_list(n+1) - y_list(n))/h; % using the data for the calculation
        
        % yp_list(n) = (f(x_list(n) + h) - f(x_list(n)))/h; % forward difference formula
        % fprintf("Value %d, is %f\n", n, yp_list(n));      % for debugging
    end
    yp_list(sizeOfLists) = NaN;
end

function yp_list = firstOrderBackward(x_list, y_list)
    sizeOfLists = length(x_list);  % record of the size of your lists
    %yp_list = zeros(sizeOfLists);  % make a list

    for n = 2:sizeOfLists
        h = x_list(n) - x_list(n-1);
        yp_list(n) = (y_list(n) - y_list(n-1))/h;
    end
    yp_list(1) = NaN;

end

function yp_list = secondOrderForward(x_list, y_list)
    sizeOfLists = length(x_list);  % record of the size of your lists
    %yp_list = zeros(sizeOfLists);  % make a list

    for n = 1:sizeOfLists-2
        h = x_list(n+1) - x_list(n);
        yp_list(n) = ((-1)*(3)*y_list(n) + 4*y_list(n+1) - y_list(n+2))/(2*h);
    end
    yp_list(sizeOfLists)   = NaN;
    yp_list(sizeOfLists-1) = NaN;
end

function yp_list = secondOrderBackward(x_list, y_list)
    sizeOfLists = length(x_list);  % record of the size of your lists
    %yp_list = zeros(sizeOfLists);  % make a list

    for n = 3:sizeOfLists
        h = x_list(n) - x_list(n-1);
        yp_list(n) = ((3)*y_list(n) - 4*y_list(n-1) + y_list(n-2))/(2*h);
    end
    yp_list(1) = NaN;
    yp_list(2) = NaN;
end

function yp_list = secondOrderCenter(x_list, y_list)
    sizeOfLists = length(x_list);  % record of the size of your lists
    %yp_list = zeros(sizeOfLists);  % make a list
    

    for n = 2:sizeOfLists-1
        h = (x_list(n+1) - x_list(n-1))/2;
        yp_list(n) = (y_list(n+1) - y_list(n-1))/(2*h);
    end
    yp_list(1) = NaN;
    yp_list(sizeOfLists) = NaN;
end

function yp_list = fourthOrderCenter(x_list, y_list)
    sizeOfLists = length(x_list);  % record of the size of your lists
    %yp_list = zeros(sizeOfLists);  % make a list

    for n = 3:sizeOfLists-2
        h = (x_list(n+1) - x_list(n-1))/2;
        yp_list(n) = (y_list(n+1) - y_list(n-1))/(2*h);
    end
    yp_list(1) = NaN;
    yp_list(2) = NaN;
    yp_list( sizeOfLists ) = NaN;
    yp_list(sizeOfLists-1) = NaN;
end