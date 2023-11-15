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
  ScaVector p(size); ScaVector r(size); ScaVector r_test(size);
  ScaVector Ar(size); ScaVector Ap(size);
  double beta = 0; double alpha_ = 0;
  
  //Gradient Conjugate solver
  int it = 0;
  ScaVector Au = A*u;
  exchangeAddInterfMPI(Au, mesh);
  r = b - Au;
  p = r;
  //r_test = r;
  
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
    u += alpha_*p;

      //Update r (--> r = r - alpha_*Ap;)
    r = r - alpha_*Ap;
      
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
    if(((it % (maxit/1000)) == 0)){
       if(myRank == 0){
        //cout << "\nAprès la boucle" <<endl;
        cout << "   [" << it << "] residual: " << norm_r/norm_r_0 << endl;
        //cout << "   [" << it << "] alpha = " << alpha_ << endl;
        //cout << "   [" << it << "] beta = " << beta << endl;
        //cout << "norm_r_0 = " << norm_r_0 << endl;
       }
    }

    it++;
  }

if(myRank == 0){
    cout << "   -> final iteration: " << it << " (prescribed max: " << maxit << ")" << endl;
    cout << "   -> final residual: " << norm_r/norm_r_0 << " (prescribed tol: " << tol << ")" << endl;
  
  }
}
