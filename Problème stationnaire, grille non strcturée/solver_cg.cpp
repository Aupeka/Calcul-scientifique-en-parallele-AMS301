#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Solution of the system Au=b with gradient conjugate
//================================================================================

void gradient_conjugate(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh, double tol, int maxit)
{
  if(myRank == 0)
    cout << "== gradient conjugate" << endl;

  // Compute the solver parameters
  int size = A.rows();
  ScaVector p(size); ScaVector r(size);
  ScaVector Ar(size); ScaVector Ap(size);
  double beta; double alpha_;
  
  //Gradient Conjugate solver
  int it = 0;
  ScaVector Au = A*u;
  exchangeAddInterfMPI(Au, mesh);
  r = b - Au;
  p = r;
  
  //Residu initialization
  double norm_r = 1e2;
  double norm_r_0 = norm_2_glo(r,mesh);

  while (norm_r/norm_r_0 > tol && it < maxit){
    //==============================================
    // 1. Updates
    //==============================================

      //Update alpha_
    Ap = A*p;
    exchangeAddInterfMPI(Ap, mesh);
    alpha_ = produit_scalaire_glo(r,p,mesh)/produit_scalaire_glo(Ap,p,mesh);

      //Update field
    u = u + alpha_*p;

      //Update r
    r = r - alpha_*Ap;
      /*Au = A*u;
      exchangeAddInterfMPI(Au, mesh);
      r = b - Au;*/   
    
      //Update beta
    Ar = A*r;
    exchangeAddInterfMPI(Ar, mesh);

    beta = -produit_scalaire_glo(Ar,p,mesh)/produit_scalaire_glo(Ap,p,mesh);

      //Update p
    p = r + beta*p;

      //Update norm_r
    norm_r = norm_2_glo(r,mesh);

    //==============================================
    // 2. Display
    //==============================================
    if(((it % (maxit/100)) == 0)){
       if(myRank == 0){
        cout << "   [" << it << "] residual: " << norm_r/norm_r_0 << endl;
        //cout << "alpha = " << alpha_ << endl;
        //cout << "beta = " << beta << endl;
        //cout << "norm_r = " << norm_r << endl;
       }
    }

    it++;
  }

if(myRank == 0){
    cout << "   -> final iteration: " << it << " (prescribed max: " << maxit << ")" << endl;
    cout << "   -> final residual: " << norm_r/norm_r_0 << " (prescribed tol: " << tol << ")" << endl;
  
  }

}
