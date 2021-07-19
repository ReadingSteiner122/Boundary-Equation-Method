function y = intg(t, xi, eta, xk, yk, nkx, nky, lk)
y = (nkx*(xk - t*lk*nky - xi) + nky*(yk + t*lk*nkx - eta))...
./((xk - t*lk*nky - xi).^2 + (yk + t*lk*nkx - eta).^2);