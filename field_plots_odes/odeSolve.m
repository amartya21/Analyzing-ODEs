% ODE solve

%% QUESTION 1
x = linspace(-4,4,1000);
y1 = x.^8;
y2 = 4.^x;
figure
plot(x,y1);
hold on
plot(x,y2);
axis([-4 3 -4 10]); %shifts the view of the 2D graph
title("Graphs: of x^8 and 4^x")
h=@(x) x.^8-4.^x;   %difference of the function
p1=fzero(h,-1)      % the first point of intersection 
p2=fzero(h,1)       % the second point of intersection

%% QUESTION 2
%(a)
syms x;     % declare symbolic variable
disp("The derivative for (a) is: ")
diff((x.^3/(x.^2 + 1)),x)

%(c)
disp("The third derivative for (c) is: ")
diff(diff(diff(atan(x))))

%(e)
disp("The derivative for (e) is: ")
diff((exp(x.*log(x))),x)

%% QUESTION 3
%(a)
syms t      %symbolic variable
f =@(t) t.^6 - 4*t.^4 - 2*t.^3 + 3*t.^2 + 2*t;  %function definition
figure
ezplot(f,[-3/2,5/2]);
title('Q 10: t^6 - 4*t^4 - 2*t^3 + 3*t^2 + 2*t')
%(b)
% 2 maximas and 2 minimas are visible
axis([-1.5 2.5 -8 5]);
disp("Q10(b): 2 maximas and 2 minimas are visible")

%(c)
derivativef=diff(f,t)
derivativef = @(t) 6*t^5 - 16*t^3 - 6*t^2 + 6*t + 2;
disp('Inflection points @ -1,0,1 and 2');
fzero(derivativef,-1)
fzero(derivativef,0)
fzero(derivativef,1)
fzero(derivativef,2)

%(d)
figure
ezplot(diff(derivativef,t),[-1.2, -0.8]);
title("Second derivative of Q10")

%% QUESTION 4
%(a)
figure
syms y x
func = 3.*y+y.^3-x.^3-5
ezplot(func, [-10,10])
title("Q13 a")

%(b)
figure
[x,y] = meshgrid(-10:0.1:10, -10:0.1:10);
z = 3.*y+y.^3-x.^3;
contour(x,y,z,[-2,2])
hold on
contour(x,y,z,[0,0])
hold on
contour(x,y,z,[2,2])
hold on
contour(x,y,z,[5,5])
hold on
contour(x,y,z,[8,8])
title("Q13 b")
hold off
%(c)
% func = 0 when x=y=1
figure
[x,y] = meshgrid(0.1:0.1:3);
z = x.*log(y)+y.*log(x);
contour(x,y,z,[0,0])
title("Q13 c")