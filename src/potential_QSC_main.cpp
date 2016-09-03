#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <cassert>

#include "potential_QSC_class.h"
#include "potential_QSC_parameters.h"

/* 
 * This is an implementation of the Quantum Sutton-Chen Potential,
 * which is an EAM type potential with the following functional form:
 * EAM functional: F_i(rho_i) = c*sqrt(rho_i)
 * EAM density:    rho_i(r_ij) = sum_(i!=j) (a/r_ij)^m
 * Pair potential: V(r_ij) = (a/r_ij)^n
 */


void QSC::QSC_constructor(void)
{
    cutoff = 10.0;
    verlet_skin = 0.5;
    init = false;
    prev_num_atoms = 0;
    vlist_updates = 0;

    int i=0;
    while (true) {
        if (qsc_default_params[i].Z == -1) break;
        qsc_params.push_back(qsc_default_params[i]);
        i++;
    }

}

void QSC::QSC_destructor(void)
{
    cleanMemory();
}

QSC::QSC()
{
    cutoff = 10.0;
    verlet_skin = 0.5;
    init = false;
    prev_num_atoms = 0;
    vlist_updates = 0;

    int i=0;
    while (true) {
        if (qsc_default_params[i].Z == -1) break;
        qsc_params.push_back(qsc_default_params[i]);
        i++;
    }

}

QSC::~QSC()
{
    cleanMemory();
}

void QSC::initialize(long N, const double *R, const int *atomicNrs,
                     const double *box)
{
        /* Allocate memory */
        natomstoclear = N; //this is ugly but i need this in cleanMemory
        rho = new double[N];
        sqrtrho = new double[N];
        vlist = new int*[N];
        nlist = new int[N];
        distances = new struct distance*[N];
        oldR = new double[3*N];
        V = new double*[N];
        phi = new double*[N];
        for (int i=0;i<N;i++) {
            vlist[i] = new int[N];
            distances[i] = new struct distance[N];
            V[i] = new double[N];
            phi[i] = new double[N];
        }

        init = true;
}


void QSC::cleanMemory()
{
    if (init==true) {
        delete rho;
        delete sqrtrho;
        delete oldR;
        delete nlist;
        for (int i=0; i<natomstoclear; i++) {
            delete vlist[i];
            delete distances[i];
            delete V[i];    
            delete phi[i];
        }
        delete vlist;
        delete distances;
        delete V;    
        delete phi;

        init = false;
    }
}

void QSC::new_vlist(long N, const double *R, const double *box)
{
    double rv = cutoff+verlet_skin;

    for (int i=0; i<N; i++) {
        nlist[i] = 0;
        for (int j=i+1;j<N; j++) {
            calc_distance(box, R, i, R, j, &distances[i][j]);
            if (distances[i][j].r <= rv) {
                vlist[i][nlist[i]] = j;
                nlist[i] += 1;
            }
        }
    }

    for (int i=0; i<3*N; i++) {
        oldR[i] = R[i];
    }

    vlist_updates++;
}

void QSC::update_distances(long N, const double *R, const double *box) 
{
    for (int i=0; i<N; i++) {
        for (int k=0; k<nlist[i]; k++) {
            int j = vlist[i][k];
            calc_distance(box, R, i, R, j, &distances[i][j]);
        }
    }
}

bool QSC::verlet_needs_update(long N, const double *R, const double *box)
{
    distance diff;
    double dist_max1=0;
    double dist_max2=0;
    double dist_sum;
    
    for (int i=0; i<N; i++) {
        calc_distance(box, oldR, i, R, i, &diff);
        if(diff.r > dist_max1) {
            dist_max2 = dist_max1;
            dist_max1 = diff.r;
        } else if(diff.r > dist_max2) {
            dist_max2 = diff.r;
        }
        dist_sum = dist_max1+dist_max2;
        if (dist_sum > verlet_skin) {
            return true;
            break;
        }
    }
    return false;
}

void QSC::energy(long N, const double *R, const int *atomicNrs, double *U,
                 const double *box) {
    *U = 0.0;
    for(int i=0;i<N;i++){
        rho[i] = 0.0;
    }

    /* Calculate the local density (rho[i]) for each atom 
     * and the potential energy (U). */
    int prev_i_atomic_number=-1;
    int prev_j_atomic_number=-1;
    for (int i=0; i<N; i++) {
        double pair_term=0.0;
        qsc_parameters p_ii;

        /* Get the parameters for element i */
        if (prev_i_atomic_number != atomicNrs[i]) {
            p_ii = get_qsc_parameters(atomicNrs[i], atomicNrs[i]);
        }
        prev_i_atomic_number = atomicNrs[i];

        for (int k=0; k<nlist[i]; k++) {
            qsc_parameters p_ij, p_jj;
            int j = vlist[i][k];
            if (distances[i][j].r > cutoff) continue;

            if (prev_j_atomic_number != atomicNrs[j]) {
                p_ij = get_qsc_parameters(atomicNrs[i], atomicNrs[j]);
                p_jj = get_qsc_parameters(atomicNrs[j], atomicNrs[j]);
            }
            prev_j_atomic_number = atomicNrs[j];

            /* Take care of density */
            phi[i][j] = pair_potential(distances[i][j].r, p_jj.a, p_jj.m);
            rho[i] += phi[i][j];
            phi[j][i] = pair_potential(distances[i][j].r, p_ii.a, p_ii.m);
            rho[j] += phi[j][i];

            /* Repulsive pair term */
            V[i][j] = p_ij.epsilon*
                      pair_potential(distances[i][j].r, p_ij.a, p_ij.n);

            pair_term += V[i][j];
        }
        sqrtrho[i] = sqrt(rho[i]);
        double embedding_term = p_ii.c * p_ii.epsilon * sqrtrho[i];
        *U += pair_term - embedding_term;
    }

}

void QSC::force(long N, const double *R, const int *atomicNrs, double *F,
                double *U, const double *box)
{
    if (init==false) {
        initialize(N, R, atomicNrs, box);
        new_vlist(N, R, box);
        prev_num_atoms = N;
    }else if (N != prev_num_atoms) {
        cleanMemory();
        initialize(N, R, atomicNrs, box);
        new_vlist(N, R, box);
        prev_num_atoms = N;
    }else if (verlet_needs_update(N, R, box)) {
        new_vlist(N, R, box);
    }else{
        update_distances(N, R, box);
    }

    energy(N, R, atomicNrs, U, box);

    /* Zero out Forces */
    for(int i=0;i<3*N;i++){
        F[i] = 0.0;
    }

    /* Forces Calculation */
    int prev_i_atomic_number=-1;
    int prev_j_atomic_number=-1;
    for (int i=0; i<N; i++) {
        qsc_parameters p_ii;
        if (prev_i_atomic_number != atomicNrs[i]) {
            p_ii = get_qsc_parameters(atomicNrs[i], atomicNrs[i]);
        }
        prev_i_atomic_number = atomicNrs[i];
        for (int k=0; k<nlist[i]; k++) {
            qsc_parameters p_ij, p_jj;
            int j = vlist[i][k];

            double r_ij = distances[i][j].r;
            if (r_ij > cutoff) continue;

            if (prev_j_atomic_number != atomicNrs[j]) {
                p_ij = get_qsc_parameters(atomicNrs[i], atomicNrs[j]);
                p_jj = get_qsc_parameters(atomicNrs[j], atomicNrs[j]);
            }
            prev_j_atomic_number = atomicNrs[j];

            double Fij;
            Fij  = p_ij.n*V[i][j];
            Fij -= p_ii.epsilon*p_ii.c*p_jj.m*0.5*(1.0/sqrtrho[i])*phi[i][j];
            Fij -= p_jj.epsilon*p_jj.c*p_ii.m*0.5*(1.0/sqrtrho[j])*phi[j][i];
            Fij /= r_ij;

            double Fijx = Fij * distances[i][j].d[0]/r_ij;
            double Fijy = Fij * distances[i][j].d[1]/r_ij;
            double Fijz = Fij * distances[i][j].d[2]/r_ij;

            F[3*i]   += Fijx;
            F[3*i+1] += Fijy;
            F[3*i+2] += Fijz;
            F[3*j]   -= Fijx;
            F[3*j+1] -= Fijy;
            F[3*j+2] -= Fijz;
        }
    }
}

double QSC::dpowi(double x, unsigned n)
{
    double p;
    double r;

    p = x;
    r = 1.0;

    while (n>0) {
        if (n%2 == 1)
            r *= p;
        p *= p;
        n /= 2;
    }

    return r;
}

double QSC::pair_potential(double r, double a, double n)
{
    double result;
    double x = a/r;

    if ( (n-floor(n)) == 0.0  && n>0) {
        result = dpowi(x, (unsigned int)n);
    }else{
        result = pow(x, n);
    }

    return result;
}

void QSC::calc_distance(const double *box, const double *R1, int i,
                        const double *R2, int j, struct distance *d)
{
    double diffRX = R1[3*i]   - R2[3*j];
    double diffRY = R1[3*i+1] - R2[3*j+1];
    double diffRZ = R1[3*i+2] - R2[3*j+2];

    /* Orthogonal PBC */
    diffRX = diffRX-box[0]*floor(diffRX/box[0]+0.5); 
    diffRY = diffRY-box[4]*floor(diffRY/box[4]+0.5);
    diffRZ = diffRZ-box[8]*floor(diffRZ/box[8]+0.5);
    
    d->r = sqrt(diffRX*diffRX+diffRY*diffRY+diffRZ*diffRZ);
    d->d[0] = diffRX;
    d->d[1] = diffRY;
    d->d[2] = diffRZ;
}

void QSC::set_verlet_skin(double dr)
{
    assert(dr > 0.0);
    verlet_skin = dr;
}

void QSC::set_cutoff(double c)
{
    assert(c > 0.0);
    cutoff = c;
}

double QSC::get_cutoff(void)
{
    return cutoff;
}


//Either adds or modifies the parameters list
void QSC::set_qsc_parameter(int Z, double n, double m, double epsilon,
                            double c, double a)
{
    qsc_parameters p;
    p.Z = Z;
    p.n = n;
    p.m = m;
    p.epsilon = epsilon;
    p.c = c;
    p.a = a;

    bool match=false;
    for (unsigned int i=0;i<qsc_params.size(); i++) {
        if (qsc_params[i].Z == p.Z) {
            qsc_params[i] = p;
            match=true;
            break;
        }
    }
    if (!match) {
        qsc_params.push_back(p);
    }
}

QSC::qsc_parameters QSC::get_qsc_parameters(int element_a, int element_b)
{
    int ia=-1, ib=-1;

    for (unsigned int i=0;i<qsc_params.size(); i++) {
        if (element_a == qsc_params[i].Z) {
            ia = i;
        }

        if (element_b == qsc_params[i].Z) {
            ib = i;
        }

        if (ia != -1 && ib != -1) break;
    }

    if (ia == -1) {
        printf("ERROR: QSC doesn't have parameters for element %i\n",
               element_a);
        throw 1;
    } else if (ib == -1) {
        printf("ERROR: QSC doesn't have parameters for element %i\n", 
               element_b);
        throw 1;
    }

    /* Mixing rules */
    if (ia==ib) {
        return qsc_params[ia];
    }



    qsc_parameters p;
    p.epsilon = sqrt(qsc_params[ia].epsilon * 
                     qsc_params[ib].epsilon);
    p.a = 0.5*(qsc_params[ia].a + qsc_params[ib].a);
    p.m = 0.5*(qsc_params[ia].m + qsc_params[ib].m);
    p.n = 0.5*(qsc_params[ia].n + qsc_params[ib].n);
    p.c = qsc_params[ia].c;
    p.Z = 0;

    return p;
}
