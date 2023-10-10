// Compilation:
//   mpicxx Jacobi_para.cpp
// Execution:
//   mpirun -np 4 ./a.out

# define M_PI = 3.14159265358979323846  /* pi */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
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

//Algorithm parameters  double L = 1e6;
    double L = 1e5;
  double epsilon = 1e-5;

// Problem parameters
  double a = 1.;
  double b = 1.;
  double alpha = 0.5;
  double U_0 = 0;
  int Nx = 50;
  int Ny = 50;

//Coefficients
  double dx = a/(Nx+1);
  double dy = b/(Ny+1);
  double coeff = (dx*dx*dy*dy)/(2*(dx*dx + dy*dy));

// Memory allocation + Initial solution + Boundary conditions
  vector<double> f((Nx+2)*(Ny+2));
  vector<double> sol((Nx+2)*(Ny+2));
  vector<double> solNew((Nx+2)*(Ny+2));

//MPI Initialization

MPI_Init(&argc, &argv);

  int nbTasks;
  int myRank;
  MPI_Comm_size(MPI_COMM_WORLD, &nbTasks);
  MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
  MPI_Barrier(MPI_COMM_WORLD); 

  // i
  int i_start = myRank * ((Nx+2)/nbTasks);
  int i_end = min((myRank+1) * ((Nx+2)/nbTasks), Nx+2);

  // j
  //int j_start = myRank * ((Ny+2)/nbTasks);
  //int j_end = min((myRank+1) * ((Ny+2)/nbTasks), (Ny+2));

  // n
  int n_start = myRank * ((Nx+2)*(Ny+2)/nbTasks);
  int n_end = min((myRank+1) * ((Nx+2)*(Ny+2)/nbTasks), (Nx+2)*(Ny+2));

//IC

  for (int n=0; n < (Nx+2)*(Ny+2); n++){
    sol[n] = U_0;
  }

//LC, of problem
  for (int j = 0; j < Ny+2; j++){
    sol[0 + (Nx+2)*j] = U_0*(1+alpha*V((Nx+2)*j*dy,b));
  }
  
  /* Validation process - begin*/
  
  //Validation for u = sin(2*pi*x)sin(2*pi*y), therefore f = 8*pi^2*sin(2*pi*x)sin(2*pi*y)
  
  //Memory allocation
  vector<double> sol_th((Nx+2)*(Ny+2));
  
  //f and u_solution
  for (int i = 0; i<Nx+2;i++){
    for (int j=0; j<Ny+2;j++){
      f[i + (Nx+2)*j] = -8*M_PI*M_PI*sin(2*M_PI*i*dx)*sin(2*M_PI*j*dy);
      sol_th[i + (Nx+2)*j] = sin(2*M_PI*i*dx)*sin(2*M_PI*j*dy);
    }
  }
  
  /* Validation process - end */


  /* Time loop - begin*/

  //Initiallization
  int l = 0;
  vector<double> res((Nx+2)*(Ny+2));
  double n_res = MAXFLOAT;
  double n_res_para = MAXFLOAT;
  double n_res_0 = MAXFLOAT;

  // Loop

    // l = 0
  while (l < 1){
    // Spatial loop
    for (int i=1; i<Nx+1;i++){
        for (int j=1; j<Ny+1; j++){
            solNew[i+(Nx+2)*j + 1] = coeff * ((sol[i+1+(Nx+2)*j+1]+sol[i-1+(Nx+2)*j + 1])/(dx*dx) + (sol[i+(Nx+2)*(j+1) + 1]+sol[i+(Nx+2)*(j-1) + 1])/(dy*dy) - f[i + (Nx+2)*j + 1]);     
        }
    }
    // Swap pointers
    sol.swap(solNew);
    
    //Residu
    for (int n = 0; n < (Nx+2)*(Ny+2); n++){
      res[n] = sol[n] - solNew[n];
    }
    n_res_0 = l2_norm(res);

    //Incrementation
    ++l;
  }

  // Check time
  double timeInit = MPI_Wtime();

  n_res_para = n_res_0;

    // l > 1
  while((l<=L)&&(n_res/n_res_0 > epsilon)){

    
  // Spatial loop
    for (int i=1 + i_start; i< i_end - 1; i++){
        for (int j=1; j<Ny+1; j++){
            solNew[i+(Nx+2)*j + 1] = coeff * ((sol[i+1+(Nx+2)*j+1]+sol[i-1+(Nx+2)*j + 1])/(dx*dx) + (sol[i+(Nx+2)*(j+1) + 1]+sol[i+(Nx+2)*(j-1) + 1])/(dy*dy) - f[i + (Nx+2)*j + 1]);     
            //res[i+(Nx+2)*j + 1] = solNew[i+(Nx+2)*j + 1] - sol[i+(Nx+2)*j + 1];
        }
    }

    // Swap pointers
    sol.swap(solNew);
   
  //Communication phase : update of unknows 
    
    // MPI Exchanges (with left side)
    MPI_Request reqSendLeft, reqRecvLeft, reqSendLeft_res, reqRecvLeft_res;
    MPI_Barrier(MPI_COMM_WORLD);
    if (myRank > 0){
      MPI_Isend(&sol[i_start + (Nx+2)*0],   1, MPI_DOUBLE, myRank-1, 0, MPI_COMM_WORLD, &reqSendLeft);
      MPI_Irecv(&sol[i_start + (Nx+2)* 0 - 1], 1, MPI_DOUBLE, myRank-1, 0, MPI_COMM_WORLD, &reqRecvLeft);

      MPI_Isend(&res[i_start + (Nx+2)*0],   1, MPI_DOUBLE, myRank-1, 0, MPI_COMM_WORLD, &reqSendLeft_res);
      MPI_Irecv(&res[i_start + (Nx+2)* 0 - 1], 1, MPI_DOUBLE, myRank-1, 0, MPI_COMM_WORLD, &reqRecvLeft_res);
    }

    // MPI Exchanges (with right side)
    MPI_Request reqSendRight, reqRecvRight, reqSendRight_res, reqRecvRight_res;
    if (myRank < nbTasks-1){
      MPI_Isend(&sol[i_end + (Nx+2)*(Ny+2)],   1, MPI_DOUBLE, myRank+1, 0, MPI_COMM_WORLD, &reqSendRight);
      MPI_Irecv(&sol[i_end + (Nx+2)*(Ny+2) + 1], 1, MPI_DOUBLE, myRank+1, 0, MPI_COMM_WORLD, &reqRecvRight);
   
      MPI_Isend(&res[i_end + (Nx+2)*(Ny+2)],   1, MPI_DOUBLE, myRank+1, 0, MPI_COMM_WORLD, &reqSendLeft_res);
      MPI_Irecv(&res[i_end + (Nx+2)*(Ny+2) + 1], 1, MPI_DOUBLE, myRank+1, 0, MPI_COMM_WORLD, &reqRecvLeft_res);
    }

    // MPI Exchanges (check everything is send/recv)
    if(myRank > 0){
      MPI_Wait(&reqSendLeft, MPI_STATUS_IGNORE);
      MPI_Wait(&reqRecvLeft, MPI_STATUS_IGNORE);

      //MPI_Wait(&reqSendLeft_res, MPI_STATUS_IGNORE);
      //MPI_Wait(&reqRecvLeft_res, MPI_STATUS_IGNORE);

    }
    if(myRank < nbTasks-1){
      MPI_Wait(&reqSendRight, MPI_STATUS_IGNORE);
      MPI_Wait(&reqRecvRight, MPI_STATUS_IGNORE);

      //MPI_Wait(&reqSendRight_res, MPI_STATUS_IGNORE);
      //MPI_Wait(&reqRecvRight_res, MPI_STATUS_IGNORE);
    }  

    //Residu
      for (int i=1 + i_start; i< i_end - 1; i++){
        for (int j=1; j<Ny+1; j++){
            //solNew[i+(Nx+2)*j + 1] = coeff * ((sol[i+1+(Nx+2)*j+1]+sol[i-1+(Nx+2)*j + 1])/(dx*dx) + (sol[i+(Nx+2)*(j+1) + 1]+sol[i+(Nx+2)*(j-1) + 1])/(dy*dy) - f[i + (Nx+2)*j + 1]);     
            res[i+(Nx+2)*j + 1] = sol[i+(Nx+2)*j + 1] - solNew[i+(Nx+2)*j + 1];
        }
    }

    /*
    for (int n = n_start; n < n_end; n++){
      res[n] = sol[n] - solNew[n];
    }*/

    //Incrementation
    ++l;
  }

  //  Norm of res
  n_res = l2_norm(res);

  // Check time
  MPI_Barrier(MPI_COMM_WORLD);
  if(myRank == 0){
    double timeEnd = MPI_Wtime();
    cout << "Runtime: " << timeEnd-timeInit << endl;
  }

  /* Time loop - end */


  /* Solution display - begin */

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
    
    for (int n = 0; n < (Nx+2)*(Ny+2); n++){
      err[n] = sol[n] - sol_th[n];
    }

    
    if (myRank == 0){

      //Print
      cout << "***** Printing the solution *****" <<endl;

      //finding out why the program stopped
      if (l > L){
        cout << "Program stopped beacause l > L" << endl;
      }

      else {
        cout << "Program stopped because the residual value =" << n_res/n_res_0 << endl;
      }

      //cout << "- Norm of the solution:   " << l2_norm(sol) << endl;
      //cout << "- Norm of the expected solution:   " << l2_norm(u_solution) << endl;

      cout << "- Absolute error: " << l2_norm(err) << endl;

      cout << "- Norm of the residual value: " << n_res/n_res_0 << endl;


      cout << "********************************" <<endl;

      //File
      ofstream file;
      file.open("Jacobi.txt");
      for (int i=0; i<=Nx+2; i++){
            for (int j=0; j<=Ny+2; j++){
              file << pos_x[i] << ";" << pos_y[j] << ";" << sol[i + Nx*j] <<endl;
        }
      }
      file.close();

    }

  /* Solution display - end */

  MPI_Finalize();

  return 0;
}
