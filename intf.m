function y = intf (t, xi, eta, xk, yk, nkx, nky, lk)
y = log ((xk - t*lk*nky - xi).^2 + (yk + t*lk*nkx - eta).^2);