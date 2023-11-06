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

double calcul_norm_residu(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh){
    
    //Calcul du résidu
    ScaVector residu = b - A*u;
    exchangeAddInterfMPI(residu, mesh);
    removeInterfMPI(residu,mesh);
    double n_residu = norm_2(residu);
    
    //MPI exchange
    double buff = n_residu;
    MPI_Allreduce (&buff, &n_residu, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);

    //return
    return n_residu;
}

void update_residu(ScaVector& residu, SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh){
    
    //Mise à jour du résidu
    residu = b - A*u;
    exchangeAddInterfMPI(residu, mesh);
    removeInterfMPI(residu,mesh);
    //cout << "norm_test = " << norm_2(residu) << endl;
}

double erreur_l2(SpMatrix& M, ScaVector& v){
    //Error
    double size = M.rows();
    ScaVector err = v.transpose()*M*v;
    double n_err = sqrt(err(0)/size);

    //MPI Allreduce
    double buff = n_err;
    MPI_Allreduce (&buff, &n_err, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);
    
    return n_err;
}

double produit_scalaire(ScaVector& u, ScaVector& v){
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i){
        accum = accum + u(i)*v(i);
    }
    return accum;

}
