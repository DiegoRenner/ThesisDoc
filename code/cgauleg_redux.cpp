inline std::pair<Eigen::RowVectorXd, Eigen::RowVectorXd>
        cgauleg_redux(int n, double eps){
    // total number of points
    int N = n*(n+1)/2-1;
    // set constant for geometric grading
    double sigma = (sqrt(2.)-1.)*(sqrt(2.)-1.);
    // initialize boundaries for first gauss-legendre rule
    double xl = sigma;
    double xr = 1.;
    // counters for mapping from gauss-legendre
    // to composite gauss-legendre
    int ii = N-1;
    int qq = n;
    // initializing points and weights
    Eigen::RowVectorXd weights_gauleg, points_gauleg;
    Eigen::RowVectorXd weights(N), points(N);
    for (int i=1; i<n; i++){
        // generate gauss-legendre rule for Interval [0,1]
        std::tie(points_gauleg, weights_gauleg) =
                gauleg(0, 1, qq, eps);
        for (int j=qq-1; j>=0; j--){
            // map to current boundaries and assign weights and points
            weights[ii] = (xr-xl)*weights_gauleg[j];
            points[ii] = (xr-xl)*points_gauleg[j]+xl;
            ii--;
        };
        qq--;
        // set new boundaries
        xr=xl;
        xl=xl*sigma;
    };
    return std::make_pair(points, weights);
}
