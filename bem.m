
%% Program to implement the Boundary Element Method
%Pre-processing
% Read boundary coodinates, values and type of boundary condition
% from input file input.dat. Takes in 4 variables per line (x, y)
% (co-ordinates), boundary type (bt) (0-Dirichlet, 1-Neumann) and boundary
% value at midpoint (bv)
Fid = fopen('input.dat','r');
indat = fscanf (Fid,'%g%g%d%g', [4, inf]);
indat = indat';
xb = indat (:,1);
yb = indat (:,2);
bt = indat (:,3);
bv = indat(:,4);
n = length (xb) - 1;
% n = number of elements
% Find array of midpoints (xm, ym) and lengths of elements (lm), and their unit
% normal vectors
for i = 1:n
    xm(i) = 0.5d0*(xb (i) + xb(i + 1));
    ym(i) = 0.5d0*(yb(i) + yb(i + 1));
    lm(i) = sqrt((xb(i + 1) - xb(i))^2d0 + (yb(i + 1) - yb(i))^2d0);
    nx(i) = (yb(i + 1) - yb(i))/lm (i);
    ny(i) = (xb(i) - xb (i + 1))/lm(i);
end
% Processing
% Find approximations for unknown boundary values by
% – constructing matrix A and vector b
% – solving the system "Ax = b" for x
for m = 1:n
    b(m) = 0d0;
    for k = 1:n
        if(k == m)
            G = 0.0;
            F = lm(k)/(2.0*pi)*(log(lm(k)/2.0) - 1.0);
            del = 1.0;
        else
            [F, G] = findfg (xm(m), ym (m), xb(k), yb(k), nx(k), ny(k), lm(k));
            del = 0.0;
        end
        if (bt(k) == 0)
            A(m, k) = -F;
            b(m) = b(m) + bv(k)*(- G + 0.5d0*del);
        else
            A(m, k) = G - 0.5d0*del;
            b(m) = b(m) + bv(k)*F;
        end
    end
end
z = A\b';
% solve system "Ax = b" and store in z
% Assign approximate boundary values accordingly
for m = 1:n
    u (m) = (1 - bt (m))*bv (m) + bt (m)*z(m);
    q(m) = (1 - bt (m))*z(m) + bt(m)*bv (m);
end
%------ Post-processing ------
% Find value(s) at required point (s)
for j = 1:99
    y(j) = 0.01*j;
    for i = 1:99
        x(i) = 0.01*i;
        s(j, i) = 0d0;
        for k = 1:n
            [F, G] = findfg (x(i), y(j), xb(k), yb(k), nx(k), ny(k), lm(k));
            s(j, i) = s(j, i) + u(k)*G - q(k)*F;
        end
    end
end
figure(1)
surface(x, y, s,'EdgeColor','none')
hold on
% graphics of elements
for i = 1:n
    xx = [xb(i),xb(i+ 1)];
    yy = [yb(i), yb(i + 1)];
    line(xx, yy,'Color','k','LineWidth',3)
end
title ('Surface Plot of Solution')
xlabel ('x')
ylabel('y')
axis([0 1 0 1])
axis equal
