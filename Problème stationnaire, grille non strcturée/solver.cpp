#include "headers.hpp"

extern int myRank;
extern int nbTasks;

//================================================================================
// Solution of the system Au=b with Jacobi
//================================================================================

void jacobi(SpMatrix& A, ScaVector& b, ScaVector& u, Mesh& mesh, double tol, int maxit)
{
  if(myRank == 0)
    cout << "== jacobi" << endl;
  
  // Compute the solver matrices : M (matrice diag) et N (A - M) 
  int size = A.rows();
  ScaVector Mdiag(size);
  SpMatrix N(size, size);
  for(int k=0; k<A.outerSize(); ++k){
    for(SpMatrix::InnerIterator it(A,k); it; ++it){
      if(it.row() == it.col())
        Mdiag(it.row()) = it.value();
      else
        N.coeffRef(it.row(), it.col()) = -it.value();
    }
  }
  exchangeAddInterfMPI(Mdiag, mesh);//Echange les coeffs pour les noeuds des deux interfaces et on calcul le résultat
  
  // Jacobi solver
  double residuNorm = 1e2;
  int it = 0;
  
  //Residu initialization
  ScaVector residu(size);
  ScaVector residu_0(size);
  residu_0 = b - A*u;
  double residuNorm_0 = norm_2(residu_0);
  
  while (residuNorm > tol && it < maxit){
    
    // Compute N*u
    ScaVector Nu = N*u;
    exchangeAddInterfMPI(Nu, mesh); //Afficher avant-après pour voir la modif
    
    // Update field
    for(int i=0; i<size; i++){
      u(i) = 1/Mdiag(i) * (Nu(i) + b(i));
    }
    // Update residual and iterator
    if((it % 100) == 0){
      residu = b - A*u;
      residuNorm = norm_2(residu)/residuNorm_0; //Tableau noeuds d'interface donc voir comment faire : retrancher --> possibilité d'avoir des noeuds entre 3 sous domaines
      if(myRank == 0)
        cout << "   [" << it << "] residual: " << residuNorm << endl;
    }
    it++;
  }
  
  if(myRank == 0){
    cout << "   -> final iteration: " << it << " (prescribed max: " << maxit << ")" << endl;
    cout << "   -> final residual: " << residuNorm << " (prescribed tol: " << tol << ")" << endl;
  }
}
