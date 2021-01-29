double rtsafe( std::function<double(double)> fct,
               std::function<Eigen::MatrixXd(double)> fct_both,
               double x1,
               double x2,
               double tol,
               bool &root_found,
               unsigned &num_iter){
    // initialize counter, function values and boundaries
    int j;
    double df,dx,dxold,f,fh,fl;
    double temp,xh,xl,rts;
    fl = fct(x1);
    fh = fct(x2);
    // sanity checks
    if ((fl > 0.0 && fh > 0.0) || (fl < 0.0 && fh < 0.0)) {
        std::cout << "Root must be bracketed in rtsafe" << std::endl;
        if (fl > fh){
            return x1;
        }else{
            return x2;
        };
    }
    if (fl == 0.0) {
        root_found = true;
        return x1;
    }
    if (fh == 0.0) {
        root_found = true;
        return x2;
    }
    // assign boundaries s.t f(xl) < 0
    if (fl < 0.0) {
        xl=x1;
        xh=x2;
    } else {
        xh=x1;
        xl=x2;
    }
    // set first guess for root, last step, and "step before last"
    rts=0.5*(x1+x2);
    dxold=fabs(x2-x1);
    dx=dxold;
    // initialize functions and derivative value
    Eigen::MatrixXd tempM = fct_both(rts);
    f = tempM(0,0);
    df = tempM(0,1);
    // loop over max number of allowed iterations
    for (j=1;j<=MAXIT;j++) {
        // check if newton step is out of range or too slow
        // if it is, bisect
        if ((((rts-xh)*df-f)*((rts-xl)*df-f) > 0.0) || (fabs(2.0*f) > fabs(dxold*df))) {
            dxold=dx;
            dx=0.5*(xh-xl);
            rts=xl+dx;
            // if change in root is negligible assume convergence
            if (xl == rts) {
                root_found = true;
                num_iter = j;
                return rts;
            }
        } else {
            // newton step accepted
            dxold=dx;
            dx=f/df;
            temp=rts;
            rts -= dx;
            // if change in root is negligible assume convergence
            if (temp == rts) {
                root_found = true;
                num_iter = j;
                return rts;
            }
        }
        // if change in root is small enough assume convergence
        if (fabs(dx) < tol) {
            root_found = true;
            num_iter = j;
            return rts;
        }
        // one new function/derivative evaluation per iteration
        tempM = fct_both(rts);
        f = tempM(0,0);
        df = tempM(0,1);
        // assign new boundary making sure root stay bracketed
        if (f < 0.0)
            xl=rts;
        else
            xh=rts;
    }
    std::cout << "Maximum number of iterations exceeded in rtsafe" << std::endl;
    // should never be reached
    return 0.0;
};
