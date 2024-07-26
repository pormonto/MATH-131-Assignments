

%%%%%%%%%%%%%%%%%%
%%% QUESTION 6 %%%
function Sum = my_sum(n)
    Sum = 0;
    for k = 1:n
        Sum = Sum + k;
    end
    %return Sum; % purposeful error for Q5
end

fprintf("%d\n", my_sum(5  ));
fprintf("%d\n", my_sum(20 ));
fprintf("%d\n\n", my_sum(100));


%%%%%%%%%%%%%%%%%%
%%% QUESTION 7 %%%
function value = myabsolutevalue(a)
    if a < 0
        value = a * -1;
    elseif a == 0
        value = 0;
        fprintf("Warning: Your input is 0. It cannot be positive or negative.\n")
    else
        value = a;
    end
end

fprintf("%f\n", myabsolutevalue(-2.3343224));
fprintf("%f\n", myabsolutevalue(9242.23));
fprintf("%f\n\n", myabsolutevalue(0));


%%%%%%%%%%%%%%%%%%
%%% QUESTION 8 %%%
function norm = vectornorm(x)
    norm = 0;
    for k = 1:numel(x)
        k_sq = x(k)^2;
        norm = k_sq + norm;
    end
    norm = sqrt(norm);
end

fprintf("%f\n", vectornorm([1,1,1]));
fprintf("%f\n", vectornorm([1/sqrt(2), 0, 1/sqrt(2)]));
fprintf("%f\n\n", vectornorm(0:0.01:1));


%%%%%%%%%%%%%%%%%%%
%%% QUESTION 10 %%%
h_values = logspace(-20,0,200); % logspace to create an array of 200 h_values
derlimlist = [numel(h_values)]; % making an empty list to put our approximated derivatives
relative_error = [numel(h_values)]; % making an empty list to put our relative errors

for item = 1:numel(h_values)
    % disp(h_values(item)); % debugging potential logic error
    h = h_values(item);
    derlimlist(item) = (sin(1.2 + h) - sin(1.2))/h; % derivative approximation calculation
    relative_error(item) = myabsolutevalue( derlimlist(item)-cos(1.2) ) / cos(1.2); % calculated the relative error and added it to the list
end

loglog(h_values, relative_error);

% DISCUSSION %
% As h gets closer to 10^-20, or in our case, 0, the relative error has a
% cusp in which the bottom of the cusp has the lowest relative error for
% the approximation of the derivative. This minima is located at
% (9.11589 x 10^-9, 8.89569 x 10^-12).
