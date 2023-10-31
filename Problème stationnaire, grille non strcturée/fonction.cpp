#include "headers.hpp"

double norm_2(ScaVector& u) { //size peut être enlevé
    double size = u.size();
    double accum = 0.;
    for (int i = 0; i < size; ++i) {
        accum += u(i)*u(i);
    }
    //MPI.Reduc
    return sqrt(accum/size);
}

double erreur_l2(SpMatrix& M, ScaVector& v){
    double size = M.rows();
    
    //Calcul de l'erreur
    ScaVector err = v.transpose()*M*v;
    
    return sqrt(err(0)/size);
}
