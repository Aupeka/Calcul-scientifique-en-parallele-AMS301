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

  // Gradient Conjugate solver
  double p = 1e2;
  double alpha_;
  int it = 0;

  //Residu initialization
  double p_0 = calcul_residu(A,b,u, mesh);

  while (residuNorm/residuNorm_0 > tol && it < maxit){
    // Update field
    alpha_ = (p.transpose()*p)/(p.transpose()*A*p);
    u += alpha_*p;
    p = calcul_residu(A,b,u, mesh);

    it++;
  }

if(myRank == 0){
    cout << "   -> final iteration: " << it << " (prescribed max: " << maxit << ")" << endl;
    cout << "   -> final residual: " << residuNorm/residuNorm_0 << " (prescribed tol: " << tol << ")" << endl;
  }

}
