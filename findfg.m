function [F,G] = findfg (xi, eta, xk, yk, nkx, nky, lk)
F = (lk/(4.0*pi))*integral(@(t)intf(t, xi, eta, xk, yk, nkx, nky, lk),0, 1);
G = (lk/(2.0*pi))*integral(@(t)intg(t, xi, eta, xk, yk, nkx, nky, lk), 0, 1);