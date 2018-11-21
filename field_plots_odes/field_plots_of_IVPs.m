% Random field plots
%% QUESTION 1
%(a)
syms t c
sol1 = dsolve('Dy + y/t = 2','y(1) = c','t')
%(b)
tvals = [0.01,0.1,1,10]
yvals = subs(sol1,'t',tvals)
cvals = [0.8, 1, 1.2]
subs(yvals,'c',cvals)

%(c)
 for j = 0.8:0.1:1.2
    ezplot(subs(sol1,'c',j), [0,2.5])
    hold on
 end
 axis tight
 xlabel t, ylabel y
 title 'Q3, with cvalues = 0.8,0.9...1.2'
%(d)
% Irrespective of what the c value is, as t tends to infinity, the function
% goes to infinity as well. If t approaches zero from the right, then the 
% function value is dependent on the c value. So if c < 1 function tends to
% negative infinity and if c > 1, the function tends to positive infinity.
% If c = 1, then it is not defined.
%% QUESTION 2
syms y t
%(a)
eq = 'Dy*(t*exp(y) - sin(y)) + exp(y) = 0'
sol1 = dsolve(eq,'t')
f = -cos(y) - t*exp(y) %written as c = f(t,y)
%(b)
figure
ezcontour(f, [-1,4,0,3])
title 'Solutions of dy/dt = (-exp(y)/(t*exp(y) - sin(y)))'
%(c)
figure
c = subs(f, [t,y],[2,1.5]);
ezplot(f-c,[-1,4,0,3])
title 'Solutions of dy/dt with initial condition (2,1.5)'
axis([0,4,0,4])
%(d)
hold on
for j = [1,1.5,3]
    f1 = @(y) eval(subs(f,t,j)-c);
    y1 = fzero(f1,2)
    % print values for t and y
    [j, double(y1)]
    plot(j,double(y1),'o')
end
hold off
title 'Points of the solution curve of the IVP'
%% QUESTION 3
%   (a)
figure
[T, Y] = meshgrid(-2:0.2:2, -2:0.2:2);
S = -T.* Y.^3;
L = sqrt(1 + S.^2);
quiver(T,Y,1./L ,S./L,0.5), axis equal tight
xlabel 't', ylabel 'y'
title 'Direction Field for dy/dt = -ty^3'; hold off
% If t increases then t < 0, y(0) > 0; y(0) < 0 approaches negative
% infinity. When t = 0, the soultion is undefined. When t > 0, y(0) > 0 and
% y(0) < 0 both tend to zero when t increseases. There is a constant
% solution when y = 0.

%   (b)
% Solving for the explicit solution of the differential equation, we get,
% y = (1/t^2)^0.5 and y = -(1/t^2)^0.5. This is an odd function and the
% solutions are symmetric about the t-axis.

%   (c)
sol1 = dsolve('Dy =-t*y^3','y(0)=1/sqrt(c)','t')
sol1 = dsolve('Dy =-t*y^3','y(0)=-1/sqrt(c)','t')
%   (d)
figure; hold on;
syms t c;
sol = dsolve('Dy = -t*(y^3)','y(0) = 1/sqrt(c)','t');
for cval = -5:1:5
    fplot(subs(sol,'c',cval),[-3,3])
end
hold off
%f(c) tends to infinity as t tends to zero (there are asymptotes)
%(e)
% 5 different types of solution curves lying above the t-axis.(Seen in 
% figure 5 & 6) In each case t is an element of R (-inf,inf).
% For the first asymptote, interval of existence is (-inf,to itself) when
% less than 0. The second is the curve that goes from zero to infinity.
% Interval is (-inf,0). It tends to 0 at -inf and goes to inf at 0. The
% third curve goes to infinity from the right at 0. Its interval is from
% (0,inf). The curve is increasing towards zero. It tends to inf at 0 and
% goes to inf at its asymptote. Fourth curve goes to inf from the right of
% an asymptote at a value greater than 0. Its interval is (0,inf). The
% curve is increasing as it approaches the asymptote from the right. The
% fifth curve has no vertical asymptotes. Its interval is from (-inf,inf).
% Solution goes to 0 at +/- inf. (Approaches 0 from left and right).
%

%(f)
figure; hold on;
syms t c;
sol = dsolve('Dy = -t*(y^3)','y(0) = 1/sqrt(c)','t');
for cval = -5:1:5
    fplot(subs(sol,'c',cval),[-3,3])
end
hold on
[T,Y] = meshgrid([-2:.2:3], [-2:.2:3]);
S = -T.*(Y.^3);
L = sqrt(1+S.^2);
quiver(T,Y,1./L,S./L,0.5)
axis([-2 2 -2 2])
hold off
%% QUESTION 4
%(a)
syms y a
disp('roots of y -')
sol = solve((a - 1)*y -y.^3)
% y = 0 is the only real root when a <= 1. If this was the case, y^2 = -constant.
% The roots would be imaginary. 
% There would be 3 distinct real roots only when a > 1.

%(b)
a=-1;
[T,Y] = meshgrid([-1:.1:1], [-1:.1:1]);
S = (a-1).*Y-Y.^3;
L = sqrt(1+S.^2);
figure
quiver(T,Y,1./L,S./L)
title('a=-1');
xlabel('t');
ylabel('y');
axis tight
% Stable in all cases
a=0;
[T,Y] = meshgrid([-1:.1:1], [-1:.1:1]);
S = (a-1).*Y-Y.^3;
L = sqrt(1+S.^2);
figure
quiver(T,Y,1./L,S./L)
title('a=0');
xlabel('t');
ylabel('y');
axis tight
% All of the graphs appear similar and all approach 0 making it stable in
% all cases

%(c)
a=1;
[T,Y] = meshgrid([-1:.1:1], [-1:.1:1]);
S = (a-1).*Y-Y.^3;
L = sqrt(1+S.^2);
figure
quiver(T,Y,1./L,S./L)
title('a=1');
xlabel('t');
ylabel('y');
axis tight
% Stable in all cases
%(d)
a=1.5;
[T,Y] = meshgrid([-2:.2:2], [-2:.2:2]);
S = (a-1).*Y-Y.^3;
L = sqrt(1+S.^2);
figure
quiver(T,Y,1./L,S./L)
title('a=1.5');
xlabel('t');
ylabel('y');
axis tight
% The stationary points appear to be at -0.7 (stable), 0 (unstable), and
% 0.7 (which is stable)

a=2;
[T,Y] = meshgrid([-2:.2:2], [-2:.2:2]);
S = (a-1).*Y-Y.^3;
L = sqrt(1+S.^2);
figure
quiver(T,Y,1./L,S./L)
title('a=2');
xlabel('t');
ylabel('y');
axis tight
% The equilibrium points appear to be at -1 (stable), 0 (unstable), and
% 1 (stable) We can see the direction fields diverge when unstable.

%(e)
% When 'a' tends to 1, the solutions after the stationary solution 0
% bifurcates because the +-(a-1)^(0.5) becomes real and non imaginary. This
% means that now the function can take one of the parts that have now 
% become real. This makes 3 stationary points and 
% this implies 0 unstable when the parts (fields) diverge