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
  ScaVector uNum(mesh.nbOfNodes);
  ScaVector uExa(mesh.nbOfNodes);
  ScaVector f(mesh.nbOfNodes);
  for(int i=0; i<mesh.nbOfNodes; ++i){
    double x = mesh.coords(i,0);
    double y = mesh.coords(i,1);
    uNum(i) = 0.;
    uExa(i) = sin(2*M_PI*x)*sin(2*M_PI*y);//cos(4*M_PI*x)*cos(M_PI*y); // a = 1 = b
    f(i) = 8*M_PI*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y);//17*M_PI*M_PI*uExa(i); // a = 1 = b
  }
  
  Problem pbm;
  double alpha = 1;
  buildProblem(pbm,mesh,alpha,f);
  
  // 4. Solve problem
  double tol = 1e-6;
  int maxit = 1e4;
  jacobi(pbm.A, pbm.b, uNum, mesh, tol, maxit);
  
  // 5. Compute error and export fields
  ScaVector uErr = uNum - uExa;
  
    //*Export fields
  exportFieldMsh(uNum, mesh, "solNum", "benchmark/solNum.msh");
  exportFieldMsh(uExa, mesh, "solRef", "benchmark/solExa.msh");
  exportFieldMsh(uErr, mesh, "solErr", "benchmark/solErr.msh");
  
    //Compute error
  double err_l2 = erreur_l2(pbm.M, uErr); //uErr = v
  if (myRank == 0){cout << "   -> Erreur L2 : " << err_l2 << endl;}
  
  // 6. Finilize MPI
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
  
  return 0;
}
