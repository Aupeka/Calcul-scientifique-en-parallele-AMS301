#include "headers.hpp"

extern int myRank;
extern int nbTasks;

double norm_2(ScaVector u) {
    double size = u.size();
    double accum = 0.;
    for (int i = 0; i < size; ++i) {
        accum += u(i)*u(i);
    }
    return sqrt(accum/size);
}

double norm_2_glo(ScaVector u, Mesh& mesh) {
    removeInterfMPI(u,mesh);
    double n_loc = 0;
    double n_glo;

    double size = u.size();
    for (int i = 0; i < size; ++i) {
        n_loc += u(i)*u(i);
    }
    
    MPI_Allreduce (&n_loc, &n_glo, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);

    return sqrt(n_glo); //size --> Plus haut
}

double produit_scalaire(ScaVector u, ScaVector& v){
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i){
        accum = accum + u(i)*v(i);
    }
    return accum;
}

double produit_scalaire_glo(ScaVector u, ScaVector v, Mesh& mesh){
    removeInterfMPI(u,mesh);
    removeInterfMPI(v,mesh);
    double ps_loc = (u.transpose()*v)(0);
    double ps_glo;

    MPI_Allreduce (&ps_loc, &ps_glo, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);

    return ps_glo; // Diviser par le size ????
}

double erreur_l2(SpMatrix& M, ScaVector v, Mesh& mesh){
    //Error
    double size = M.rows();
    removeInterfMPI(v,mesh);
    ScaVector err = v.transpose()*M*v;
    double n_err = err(0);

    //MPI Allreduce
    double buff = n_err;
    MPI_Allreduce (&buff, &n_err, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
    
    return sqrt(n_err);
}
