// integrating by tensor product quadrature rule
complex_t integral = complex_t(0.,0.);
// loop over quadrature rule for each dimension
for (unsigned int k = 0; k < N; ++k) {
	for (unsigned int l = 0; l < N; ++l) {
		// generate points and weights from existing QR
		double s = GaussQR.x(l)*(1.-GaussQR.x(k));
		double t = GaussQR.x(l);
		double w = GaussQR.x(l)*GaussQR.w(k)*GaussQR.w(l);
		// evaluate integrand symetrically
		integral += w*integrand(s,t);
		integral += w*integrand(t,s);
	}
}
