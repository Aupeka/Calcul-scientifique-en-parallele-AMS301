#include "headers.hpp"

int myRank;
int nbTasks;

int main(int argc, char* argv[])
{
  
  // 1. Initialize MPI
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);
 
  // 2. Read the mesh, and build lists of nodes for MPI exchanges, local numbering
  Mesh mesh;
  readMsh(mesh, "benchmark/mesh.msh");
  buildListsNodesMPI(mesh);
  
  // 3. Build problem (vectors and matrices)
  double alpha = 1;
  ScaVector uNum(mesh.nbOfNodes);
  ScaVector uExa(mesh.nbOfNodes);
  ScaVector f(mesh.nbOfNodes);
  for(int i=0; i<mesh.nbOfNodes; ++i){
    double x = mesh.coords(i,0);
    double y = mesh.coords(i,1);
    uNum(i) = 0.;
    uExa(i) = cos(4*M_PI*x)*cos(M_PI*y); // a = 1 = b
    f(i) = (17*M_PI*M_PI+alpha)*uExa(i); // /!\ Pbm de la chaleur
  }
  
  Problem pbm;
  buildProblem(pbm,mesh,alpha,f);
  
  // 4. Solve problem
  double tol = 1e-6;
  int maxit = 1e5;
  //jacobi(pbm.A, pbm.b, uNum, mesh, tol, maxit);
  gradient_conjugate(pbm.A, pbm.b, uNum, mesh, tol, maxit);
  //gradient_conjugate_seq(pbm.A, pbm.b, uNum, mesh, tol, maxit);
  
  // 5. Compute error and export fields
  ScaVector uErr = uNum - uExa;
  
    //*Export fields
  exportFieldMsh(uNum, mesh, "solNum_u", "benchmark/solNum.msh");
  exportFieldMsh(uExa, mesh, "solRef", "benchmark/solExa.msh");
  exportFieldMsh(uErr, mesh, "solErr", "benchmark/solErr.msh");
  
    //Compute error
  double err_l2 = erreur_l2(pbm.M, uErr, mesh); //uErr = v
  if (myRank == 0){cout << "   -> Erreur L2 : " << err_l2 << endl;}
  
  // 6. Finilize MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}
