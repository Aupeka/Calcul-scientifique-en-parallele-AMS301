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
  int it = 0;
  ScaVector Nu;
  ScaVector Au = A*u;
  exchangeAddInterfMPI(Au, mesh);
  
  //Residu initialization
  ScaVector residu = b - Au;
  double residuNorm = 1e2;
  double residuNorm_0 = norm_2_glo(residu,mesh);
  
  // Check time
  double timeInit = MPI_Wtime();

  if(((it % (maxit/10)) == 0)){
       //if(myRank == 0)
        cout << "   [" << it << "] test: " << A << endl;
    }

  while (residuNorm/residuNorm_0 > tol && it < maxit){
    
    // Compute N*u
    Nu = N*u;
    exchangeAddInterfMPI(Nu, mesh); //Afficher avant-après pour voir la modif
    
    // Update field
    for(int i=0; i<size; i++){
      u(i) = 1/Mdiag(i) * (Nu(i) + b(i));
    }

    // Update residual and iterator
    Au = A*u;
    exchangeAddInterfMPI(Au, mesh);

    residu = b - Au;
    residuNorm = norm_2_glo(residu, mesh);

    if(((it % (maxit/10)) == 0)){
       if(myRank == 0)
        cout << "   [" << it << "] residual: " << residuNorm/residuNorm_0 << endl;
    }
    it++;
  }

  // Check time
  MPI_Barrier(MPI_COMM_WORLD);
  if(myRank == 0){
    double timeEnd = MPI_Wtime();
    cout << "   -> Runtime: " << timeEnd-timeInit << " s" << endl;
  }
  
  if(myRank == 0){
    cout << "   -> final iteration: " << it << " (prescribed max: " << maxit << ")" << endl;
    cout << "   -> final residual: " << residuNorm/residuNorm_0 << " (prescribed tol: " << tol << ")" << endl;
  }
}
