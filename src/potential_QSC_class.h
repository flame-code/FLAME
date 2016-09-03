#include <vector>

#ifndef QSC_STANDALONE
#include "../../Potential.h"
#endif

class QSC
#ifndef QSC_STANDALONE
: public Potential
#endif
{    
    public:
        QSC(void);
        ~QSC(void);

        void initialize() {};
        void initialize(long N, const double *R, const int *atomicNrs,
                        const double *box);
        void cleanMemory();
        void force(long N, const double *R, const int *atomicNrs,
                   double *F, double *U, const double *box);
        void set_verlet_skin(double dr);
        void set_cutoff(double c);
        void QSC_constructor(void);
        void QSC_destructor(void);
        double get_cutoff(void);
        void set_qsc_parameter(int Z, double n, double m, double epsilon,
                               double c, double a);

        long vlist_updates;

    private:
        struct qsc_parameters {
            int Z;
            double n;
            double m;
            double epsilon;
            double c;
            double a;
        };
        bool init;
        long prev_num_atoms;

        struct distance {
            double d[3];
            double r;
        };
        distance **distances;

        int natomstoclear;
        int **vlist;
        int  *nlist;
        double cutoff;
        double verlet_skin;
        double *oldR;
        double *rho;
        double *sqrtrho;
        double **V;
        double **phi;

        void energy(long N, const double *R, const int *atomicNrs,
                    double *U, const double *box);

        static const qsc_parameters qsc_default_params[];
        std::vector<qsc_parameters> qsc_params;
        
        qsc_parameters get_qsc_parameters(int a, int b);
        double dpowi(double x, unsigned n);
        double pair_potential(double r, double a, double n);
        bool verlet_needs_update(long N, const double *R, const double *box);
        void new_vlist(long N, const double *R, const double *box);
        void update_distances(long N, const double *R, const double *box);
        void calc_distance(const double *box, const double *R1, int i, 
                           const double *R2, int j, struct distance *d);
};
