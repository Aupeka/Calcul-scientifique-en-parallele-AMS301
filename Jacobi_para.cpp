// Compilation:
//   mpicxx Jacobi_para.cpp
// Execution:
//   mpirun -np 4 ./a.out

//# define M_PI = 3.14159265358979323846  /* pi */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <chrono>
#include <cstdlib>
using namespace std;

#include <mpi.h>

double l2_norm(vector<double> const& u) {
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i) {
        accum += u[i]*u[i];
    }
    return sqrt(accum/u.size());
}

double V(double y, double b){
    return 1 - cos(2*M_PI*y/b);
}

int main(int argc, char* argv[]){

//Algorithm parameters
  double L = 1e4;
  double epsilon = 1e-2;

// Problem parameters
  double a = 1.;
  double b = 1.;
  double alpha = 0.5;
  double U_0 = 0;
  int Nx = 10;
  int Ny = 10;

//Coefficients
  double dx = a/(Nx+1);
  double dy = b/(Ny+1);
  double coeff = (dx*dx*dy*dy)/(2*(dx*dx + dy*dy));

// Memory allocation + Initial solution + Boundary conditions
  vector<double> f((Nx+2)*(Ny+2));
  vector<double> sol((Nx+2)*(Ny+2));
  vector<double> solNew((Nx+2)*(Ny+2));

//IC

  for (int i=0; i < (Nx+2)*(Ny+2); i++){
    sol[i] = U_0;
    f[i] = 0.; //Nulle dans notre cas
  }

  
//LC, of problem
  for (int j = 0; j < Ny+2; j++){
    sol[0 + (Nx+2)*j] = U_0*(1+alpha*V((Nx+2)*j*dy,b));
  }
  
  /* ************* Validation process - begin ************* */
  //Validation for u = sin(2*pi*x)sin(2*pi*y)
  
  //Memory allocation
  vector<double> sol_th((Nx+2)*(Ny+2));
  
  //f and u_solution
  for (int i = 0; i<Nx+2;i++){
    for (int j=0; j<Ny+2;j++){
      f[i + (Nx+2)*j] = -8*M_PI*M_PI*sin(2*M_PI*i*dx)*sin(2*M_PI*j*dy);
      sol_th[i + (Nx+2)*j] = sin(2*M_PI*i*dx)*sin(2*M_PI*j*dy);
    }
  }
  /* ************* Validation process - end ************* */


  /* ************* Time loop - begin ************* */

  //Initiallization
  int l = 0;
  double n_res = 1e4;
  double n_res_0 = 0;

  // Loop

    // l = 0
  while (l < 1){
    // Spatial loop
    for (int i=1; i<=Nx; i++){
        for (int j=1; j<=Ny; j++){
            solNew[i+(Nx+2)*j + 1] = coeff * ((sol[i+1+(Nx+2)*j+1]+sol[i-1+(Nx+2)*j + 1])/(dx*dx) + (sol[i+(Nx+2)*(j+1) + 1]+sol[i+(Nx+2)*(j-1) + 1])/(dy*dy) - f[i + (Nx+2)*j + 1]);     
            n_res_0 += (solNew[i+(Nx+2)*j + 1] - sol[i+(Nx+2)*j + 1])*(solNew[i+(Nx+2)*j + 1] - sol[i+(Nx+2)*j + 1]);
        }
    }
    // Swap pointers
    sol.swap(solNew);
    
    //Residu
    n_res_0 = n_res_0/sol.size();

    //Incrementation
    l++;
  }

  //MPI Initialization
  MPI_Init(&argc, &argv);
  int nbTasks;
  int myRank;
  MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Barrier(MPI_COMM_WORLD); 

  // i
  int i_start = myRank * ((Nx+2)/nbTasks);
  int i_end = min((myRank+1) * ((Nx+2)/nbTasks), Nx+2) -1;

  //buffer for residu
  double buff;

  // timeInit
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

  // l > 1
  while((l<=L)&&(n_res/n_res_0 > epsilon)){

    /* ************* Communication phase : update of unknows ************* */
    
    // MPI Exchanges (with left side)
    MPI_Request reqSendLeft, reqRecvLeft;
    //MPI_Barrier(MPI_COMM_WORLD);
    if (myRank > 0){
      MPI_Isend(&sol[i_start*(Ny+2)],   Ny+2, MPI_DOUBLE, myRank-1, 0, MPI_COMM_WORLD, &reqSendLeft);
      MPI_Irecv(&sol[(i_start-1)*(Ny+2)], Ny+2, MPI_DOUBLE, myRank-1, 0, MPI_COMM_WORLD, &reqRecvLeft);
    }

    // MPI Exchanges (with right side)
    MPI_Request reqSendRight, reqRecvRight;
    if (myRank < nbTasks-1){
      MPI_Isend(&sol[i_end*(Ny+2)],   Ny+2, MPI_DOUBLE, myRank+1, 0, MPI_COMM_WORLD, &reqSendRight);
      MPI_Irecv(&sol[(i_end+1)*(Ny+2)], Ny+2, MPI_DOUBLE, myRank+1, 0, MPI_COMM_WORLD, &reqRecvRight);
    }

    // MPI Exchanges (check everything is send/recv)
    if(myRank > 0){
      MPI_Wait(&reqSendLeft, MPI_STATUS_IGNORE);
      MPI_Wait(&reqRecvLeft, MPI_STATUS_IGNORE);
    }
    if(myRank < nbTasks-1){
      MPI_Wait(&reqSendRight, MPI_STATUS_IGNORE);
      MPI_Wait(&reqRecvRight, MPI_STATUS_IGNORE);
    }

    /* ************* Communication phase : update of unknows ************* */
    
    //n_res re-init
    n_res = 0;

    // Spatial loop
    for (int i=i_start; i<=i_end; ++i){
        for (int j=1; j<=Ny; j++){
            //solNew[i+(Nx+2)*j + 1] = coeff * ((sol[i+1+(Nx+2)*j+1]+sol[i-1+(Nx+2)*j + 1])/(dx*dx) + (sol[i+(Nx+2)*(j+1) + 1]+sol[i+(Nx+2)*(j-1) + 1])/(dy*dy) - f[i + (Nx+2)*j + 1]);     
            solNew[i*(Ny+2) + j] = (-(sol[(i+1)*(Ny+2) + j] + sol[(i-1)*(Ny+2) + j])/(dx*dx) - (sol[i*(Ny+2) + j+1] + sol[i*(Ny+2) + j-1])/(dy*dy) + f[i*(Ny+2) + j])/(-2/(dx*dx)-2/(dy*dy));
            n_res += (solNew[i*(Ny+2) + j] - sol[i*(Ny+2) + j])*(solNew[i*(Ny+2) + j] - sol[i*(Ny+2) + j]);
        }
    }

    // Swap pointers
    sol.swap(solNew);
    
    //Residu
    n_res = n_res/sol.size();
    buff = n_res;
    MPI_Allreduce (&buff, &n_res, 1, MPI_DOUBLE , MPI_SUM, MPI_COMM_WORLD); //n_res has the right value
    
    //Incrementation
    ++l;
  }

  // timeEnd
  std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

  /* ************* Time loop - end ************* */


  /* ************* Solution display - begin ************* */

  //Def of pos
  vector<double> pos_x((Nx+2));
  vector<double> pos_y((Ny+2));

  for (int i = 0; i<Nx+2; ++i){
    pos_x[i] = i*dx;
  }

  for (int j=0; j<Ny+2; ++j){
    pos_y[j] = j*dy;
  }

  // Show results
  
    //Absolute error 
    vector<double> err((Nx+2)*(Ny+2));
    
    for (int i = 0; i < (Nx+2)*(Ny+2); i++){
      err[i] = sol[i] - sol_th[i];
    }

    
    if (myRank==0){
    
      //Print
      cout << "***** Printing Jacobi's solution *****" <<endl;

      cout << "Runtime: " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())/100 << "[ms]"<< endl;

      //finding out why the program stopped
      if (l > L){
        cout << "Program stopped beacause l = " << l << "> L" << endl;
      }

      else {
        cout << "Program stopped because the residual value = " << n_res/n_res_0 << " < "<< epsilon << " (l = " << l << ")" << endl;
      }

      cout << "- Absolute error: " << l2_norm(err) << " for h = " << 1/sqrt(Nx*Ny) << endl;

      cout << "- Norm of the residual value: " << n_res/n_res_0 << endl;

      cout << "********************************" <<endl;
    }
    

   /*
  if (myRank==0){
    ofstream file;
    file.open("Jacobi_para_cholesky.dat");
    file << "***** Printing Jacobi's solution *****" <<endl;

    file << "Runtime: " << (std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count())/100 << "[ms]"<< endl;

    //finding out why the program stopped
    if (l > L){
      file << "Program stopped beacause l = " << l << "> L" << endl;
    }

    else {
      file << "Program stopped because the residual value = " << n_res/n_res_0 << "< "<< epsilon << " where res_0 = " << n_res_0 << "and l = " << l << endl;
    }

    file << "- Absolute error: " << l2_norm(err) << " for h = " << 1/sqrt(Nx*Ny) << endl;

    file << "- Norm of the residual value: " << n_res/n_res_0 << endl;

    file << "********************************" <<endl;
    }
    
    //File
    ofstream file;
    file.open("Jacobi.dat");
    for (int i=0; i<=Nx+2; i++){
          for (int j=0; j<=Ny+2; j++){
            file << pos_x[i] << " " << pos_y[j] << " " << sol[i + (Nx+2)*j] <<endl;
      }
      file << endl;
    }
    file.close();
  */

  /* ************* Solution display - end ************* */
  
  MPI_Finalize();
  return 0;
}
