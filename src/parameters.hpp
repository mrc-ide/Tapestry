struct Parameters
{
    int K;
    double e_0, e_1;
    double v;
    double rho;

    Parameters(int K);
    Parameters(int K, double e_0, double e_1, double v, double rho);
    
    void print();
};