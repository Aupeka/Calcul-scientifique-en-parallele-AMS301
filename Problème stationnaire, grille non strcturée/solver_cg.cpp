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
  ScaVector p(size);
  ScaVector alpha_(size);
  
  //Gradient Conjugate solver
  int it = 0;
  p = b - A*u;
  exchangeAddInterfMPI(p, mesh);

  //Residu initialization
  double norm_p_0 = calcul_norm_residu(A,b,u, mesh);
  double norm_p = 1e2;
  double buff;

  while (norm_p/norm_p_0 > tol && it < maxit){
    
    // Compute alpha_
    ScaVector Ap = A*p;
    exchangeAddInterfMPI(Ap, mesh);

    for (int i = 0; i < size; i++){
        alpha_(i) = p(i)*p(i)/(p(i)*Ap(i));
    }

    // Update field
    for(int i=0; i<size; i++){
      u(i) = u(i) + alpha_(i)*p(i); // /!\ Pbm ici, des NAN sortent mais je ne sais pas pourquoi ...
    }
    cout << "norm_test = " << u << endl;

    //Update residu
    update_residu(p,A,b,u, mesh);

    //Update of norm_p
    norm_p = norm_2(p);
    buff = norm_p;
    MPI_Allreduce (&buff, &norm_p, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD);

    it++;
  }
  
if(myRank == 0){
    cout << "   -> final iteration: " << it << " (prescribed max: " << maxit << ")" << endl;
    cout << "   -> final residual: " << norm_p/norm_p_0 << " (prescribed tol: " << tol << ")" << endl;
  }

}
