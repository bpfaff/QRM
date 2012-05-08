#define RANDIN GetRNGstate();
#define RANDOUT PutRNGstate();

void frank(int *n, double *theta, double *output);
void rgig(int *n, double *r, double *s, double *p, double *k1, double *k2, double *lambda, double *chi, double *psi, double *s1, double *s2, double *xsim); 
double ef(double x, double lambda, double chi, double psi);

void SEprocExciteFunc(int *n, double *times, int *nmarks, double *marktimes, double *marks, double *beta, int *model, double *result);
double contribH(double s, double y, double gammaval, double deltaval);
double contribE(double s, double y, double gammaval, double rhoval, double deltaval);


