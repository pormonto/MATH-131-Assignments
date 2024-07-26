clc; clear;
%%%%%%
% Anonymous Functions
f  = @(x)(2*x^3 - 3*x - 1)*cos(x) - x;
fg = @(x)(2*x.^3 - 3.*x - 1).*cos(x) - x;
y  = @(x)x;
g  = @(x) exp(-x);
gp = @(x)-1*exp(-x);
h  = @(x)(atan(x)-1);
hp = @(x)(1/(x^2+1));

%%%%%%%%%%%%%%
% QUESTION 4 %
BI = bisection_method(f,-1,1,10^(-1*5),1000);
figure;
hspace = -1:0.01:1;
plot(hspace, fg(hspace), 'b');
ROOT_MATLAB = fzero(f, -1);
fprintf("The exact root according to MATLAB is %d\n", ROOT_MATLAB);
fprintf("The difference between our root and MATLAB's is %d\n\n", ROOT_MATLAB-BI(1));

% DISCUSSION %
%%%%%%%%%%%%%%
%
% i.)   The number of iterations used is 18.
%
% ii.)  The error bound is 8x10^-6.
%
% iii.) The expected rate of convergence for the bisection method O(1/2^n).
%
% iv.)  Code is above.
%
% v.)   Code is also above.
%
% NOTE: Code for the Bisection method is below.


%%%%%%%%%%%%%%
% QUESTION 5 %
FI = fixed_point_method(g,1,10^(-1*10),1000);
% plot the function g and g=x
figure;
hspace = -1:0.01:1;
plot(hspace, g(hspace), 'b'); hold on;
plot(hspace, y(hspace), 'b');


%%%%%%%%%%%%%%
% QUESTION 6 %
NE1 = Newtons_method(h,hp, 2,10^(-1*8),1000);
NE2 = Newtons_method(h,hp,-2,10^(-1*8),1000);

% DISCUSSION %
%%%%%%%%%%%%%%
%
% i.)  The first Newton's method converged, and it did so in 5 iterations.
%
% ii.) With a starting point of -2, the formula in Newton's method
% effictively finds new points that are farther and farther away from the
% root. This results in us never finding what the root is. The graphical
% way to imagine this would be the tangent line of the function at the
% given points hitting y=0 farther away at each iteration.
%
% NOTE: Code for the Newtons method is below.


%%%%%%%%%%%%%%
% QUESTION 8 %
SE = Secant_method(h,-2,2,10^(-1*8),1000);
figure;
hspace = -1:0.01:2;
plot(hspace, h(hspace), 'b');

% DISCUSSION %
%%%%%%%%%%%%%%
%
% i.)  Instead of relying on the derivative of the function at point -2,
% secant method allows us to utilize Mean Value Theorem to get a derivative
% that was between -2 and 2. In the first iteration, this meant getting
% closer to the root rather than farther away from it. In general, I
% believe that Secant method's ability to anchor itself with another point
% makes it better here.
%
% NOTE: Code for the Secant method is below.



%%%%%%%%%%%%%
% BISECTION %
function [midpoint,n,err] = bisection_method(f,a,b,tol,N)
    funcA = f(a);

    for n = 1:N
        midpoint = (a+b)/2;
        funcMid = f(midpoint);
        err = (b-a)/2;

        if abs(err) < tol  
            % check if the upper bound of the error is less than our
            % tolerance (check how this upper bound is explained it the book)
            fprintf("Root found\n");
            break

        elseif funcMid * funcA > 0  
            % its positive, so the true root in between the Midpoint and B 
            
            a = midpoint;
            funcA = f(a);
        else 
            % funcB and funcMid is negative, so the true root in between A
            % and the Midpoint

            b = midpoint;
            % funcMid is already calculated in the next round, so it is not
            % added here
        end
    end
    
    fprintf("BISECTION:\n");
    %fprintf("The vals are %f, %d, %f\n", midpoint, n, err);
    fprintf("The midpoint is %f\n", midpoint);
    fprintf("The number of iterations is %d\n", n);
    fprintf("The upper bound of err is %f\n\n", err);

end


%%%%%%%%%%%%%%%%%%%%%%%%%
% FIXED POINT ITERATION %
function [root,n,err] = fixed_point_method(g,x0,tol,N)
    for n = 1:N
        root = g(x0);   % the main fixed point calculation
        err = root-x0;  % the errpr

        if abs(err) < tol     % error check
            fprintf("Root found\n");
            break
        end
        
        x0 = root;      % replacing the result as the new input
    end

    fprintf("FIXED POINT ITERATION: \n");
    %fprintf("The vals are %f, %d, %f\n", root, n, err);
    fprintf("The midpoint is %f\n", root);
    fprintf("The number of iterations is %d\n", n);
    fprintf("The upper bound of err is %f\n\n", err);

end


%%%%%%%%%%%%%%%%%%
% NEWTONS METHOD %
function [c,n,err] = Newtons_method(f,fp,x0,tol,N)
    for n = 1:N
        
        c = x0 - (f(x0)/fp(x0)); % the Newton's method formula
        err = abs(c - x0);       % error quantity
        
        if err < tol             % error check
            fprintf("Root found\n"); 
            break
        end
        
        x0 = c;
    end

    fprintf("NEWTONS METHOD:\n");
    %fprintf("The vals are %f, %d, %f\n", c, n, err);
    fprintf("The root is %f\n", c);
    fprintf("The number of iterations is %d\n", n);
    fprintf("The upper bound of err is %f\n\n", err);

end


%%%%%%%%%%%%%%%%%
% SECANT METHOD %
function [c,n,err] = Secant_method(f,x0,x1,tol,N)
 
    for n = 2:N
        c = x1 - (f(x1)*(x1-x0))/(f(x1)-f(x0)); % The formula for Secant method
        err = abs(c - x1);                      % absolute error
        
        if err < tol                            % error check
            fprintf("Root found\n"); 
            break
        end
        
        x0 = x1;
        x1 = c;
        % pass on the values for the next generation 
        % we don't need to same the function values as shown in the
        % psuedocode because they will be updated during the calculation of
        % using the Secant method formula
        
    end
    
    fprintf("SECANT METHOD:\n");
    %fprintf("The vals are %f, %d, %f\n", c, n, err);
    fprintf("The root is %f\n", c);
    fprintf("The number of iterations is %d\n", n);
    fprintf("The upper bound of err is %f\n\n", err);

end