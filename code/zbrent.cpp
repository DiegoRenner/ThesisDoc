double zbrent( const function<double(double)> f,
               double x1,
               double x2,
               double tol,
               bool &root_found,
               unsigned &num_iter){
    // counter for iterations
    int iter;
    // initialize function values and boundaries
    double a=x1, b=x2, c=x2, d, e, min1, min2;
    double fa=f(a), fb=f(b), fc, p, q, r, s, tol1, xm;
    // sanity checks
    if ((fa > 0.0 && fb > 0.0) || (fa < 0.0 && fb < 0.0)) {
        cout << "Root must be bracketed in zbrent" << endl;
        return 0.0;
    }
    fc=fb;
    for (iter=1; iter <= MAXIT; iter++) {
        // reorient boundary for next interpolation
        if ((fb > 0.0 && fc > 0.0) || (fb < 0.0 && fc < 0.0)) {
            c = a;
            fc=fa;
            e=d=b-a;
        }
        if (fabs(fc) < fabs(fb)) {
            a=b;
            b=c;
            c=a;
            fa=fb;
            fb=fc;
            fc=fa;
        }
        // check if converged
        tol1=2.0*EPS*fabs(b)+0.5*tol;
        xm=0.5*(c-b);
        if (fabs(xm) <= tol1 || fb == 0.0){
            root_found = true;
            num_iter = iter;
            return b;
        }
        // try quadratic interpolation if
        // bounds are decreasing
        if (fabs(e) >= tol1 && fabs(fa) > fabs(fb)) {
            // compute parameters for quadratic interpolation
            s=fb/fa;
            if (a == c) {
                p=2.0*xm*s;
                q=1.0-s;
            } else {
                q=fa/fc;
                r=fb/fc;
                p=s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                q=(q-1.0)*(r-1.0)*(s-1.0);
            }
            // check if quadratic interpolation
            // would fall within bounds
            if (p > 0.0) q = -q;
            p=fabs(p);
            min1=3.0*xm*q-fabs(tol1*q);
            min2=fabs(e*q);
            if(2.0*p < (min1 < min2 ? min1 : min2)) {
                // interpolation accepted
                e=d;
                d=p/q;
            } else {
                // interpolation rejected, bisect instead
                d=xm;
                e=d;
            }
        } else {
            // convergence too slow, bounds not collapsing
            // fast enough, bisect
            d = xm;
            e = d;
        }
        // store previous best guess before
        // computing new best guess
        a=b;
        fa=fb;
        // compute new best guess
        // if step taken is too small, take
        // minimum accepted step towards 0
        if (fabs(d) > tol1)
            b += d;
        else
            b += (xm >= 0) ? fabs(tol1) : -fabs(tol1);
        // one new function evaluation per iteration
        fb = f(b);
    }
    cout << "Maximum number of iterations exceeded in zbrent" << endl;
    // should never be reached
    return 0.0;
}
