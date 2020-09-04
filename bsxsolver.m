function S = bsxsolver(deg2, N, U, LAM, alpha, h, m, n)
%#codegen
coder.inline('never')

    L = bsxfun(@times, deg2, (bsxfun(@times, N, U)));
    R = bsxfun(@times, bsxfun(@times, U', N'), deg2');
    iLAM = pinv(LAM);
    temp = pinv(iLAM-alpha*R*L);
    s = (1-alpha)*(h+alpha*(L*(temp*(R*h))));
    S = reshape(s,m,n);
end
