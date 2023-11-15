#include "headers.hpp"

extern int myRank;
extern int nbTasks;

double norm_2(ScaVector& u) {
    double size = u.size();
    double accum = 0.;
    for (int i = 0; i < size; ++i) {
        accum += u(i)*u(i);
    }
    return sqrt(accum/size);
}

double norm_2_glo(ScaVector& u, Mesh& mesh) {
    removeInterfMPI(u,mesh);
    double n_loc = (u.transpose()*u)(0);
    double n_glo;
    
    MPI_Allreduce (&n_loc, &n_glo, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);

    return sqrt(n_glo); //size --> Plus haut
}

double produit_scalaire(ScaVector& u, ScaVector& v){
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i){
        accum = accum + u(i)*v(i);
    }
    return accum;
}

double produit_scalaire_glo(ScaVector& u, ScaVector& v, Mesh& mesh){
    removeInterfMPI(u,mesh);
    removeInterfMPI(v,mesh);
    double ps_loc = (u.transpose()*v)(0);
    double ps_glo;

    MPI_Allreduce (&ps_loc, &ps_glo, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);

    return ps_glo; // Diviser par le size ????
}




double calcul_norm_residu(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh){
    
    //Calcul du résidu
    ScaVector residu = b - A*u;
    exchangeAddInterfMPI(residu, mesh);
    
    double n_residu = norm_2_glo(residu, mesh);

    /*
    //Calcul résidu
    removeInterfMPI(residu,mesh);
    double n_residu = norm_2(residu);

        //MPI exchange
    double buff = n_residu*n_residu; //Problème de Somme des racines
    MPI_Allreduce (&buff, &n_residu, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
    */

    //return
    return n_residu;
}

/*
void update_residu(ScaVector& residu, SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh){
    
    //Mise à jour du résidu
    residu = b - A*u;
    exchangeAddInterfMPI(residu, mesh);
    removeInterfMPI(residu,mesh);
}
*/

double erreur_l2(SpMatrix& M, ScaVector& v, Mesh& mesh){
    //Error
    double size = M.rows();
    removeInterfMPI(v,mesh);
    ScaVector err = v.transpose()*M*v;
    double n_err = err(0)/size;

    //MPI Allreduce
    double buff = n_err;
    MPI_Allreduce (&buff, &n_err, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
    
    return sqrt(n_err);
}
