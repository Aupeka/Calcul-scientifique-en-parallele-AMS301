// Compilation:
//   g++ projet_1.cpp
// Execution:
//   ./a.out

# define M_PI = 3.14159265358979323846  /* pi */

#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
using namespace std;

double l2_norm(vector<double> const& u) {
    double accum = 0.;
    for (int i = 0; i < u.size(); ++i) {
        accum += u[i] * u[i];
    }
    return sqrt(accum);
}

double V(double y, double b){
    return 1 - cos(2*M_PI*y/b);
}

int main(int argc, char* argv[]){

//Algorithm parameters
  double L = 500;
  double epsilon = 1e-2;

// Problem parameters
  double a = 1.;
  double b = 1.;
  double alpha = 0.5;
  double U_0 = 1;
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

//IC
  for (int i=0; i < (Nx+2)*(Ny+2); i++){
    sol[i] = U_0;
  }

  /*
  //LC, of problem
  for (int j = 0; j < Ny+2; j++){
    sol[0 + (Nx+2)*j] = U_0*(1+alpha*V((Nx+2)*j,b));
  }
  */


  /* Validation process - begin*/

  //Memory allocation
  vector<double> u_solution((Nx+2)*(Ny+2));
  
  //f and u_solution
  for (int i = 0; i<Nx+2;i++){
    for (int j=0; j<Ny+2;j++){
      f[i + (Nx+2)*j] = -2*M_PI*M_PI*cos(M_PI*i*dx)*sin(M_PI*j*dy);
      u_solution[i + (Nx+2)*j] = cos(M_PI*i*dx)*sin(M_PI*j*dy);
    }
  }

  // LC, for validation
  for (int j = 0; j < Ny+2; j++){
    sol[0 + (Nx+1)*j] = u_solution[0 + (Nx+1)*j];
    sol[Nx+2 + (Nx+1)*j] = u_solution[Nx+2 + (Nx+1)*j];
  }
  for (int i = 0; i <Nx+2; i++){
    sol[i] = u_solution[i];
    sol[i + (Nx+2)*(Ny+1)] = u_solution[i + (Nx+2)*(Ny+1)];
  }

  /* Validation process - end */


  /* Time loop - begin*/

  //Initiallization
  int l = 0;
  double res = MAXFLOAT;

  //Loop
  while((l<=L)&&(res > epsilon)){
    
    // Spatial loop
    for (int i=1; i<=Nx; i++){
        for (int j=1; j<=Ny; j++){
            solNew[i+(Nx+2)*j + 1] = coeff * ((sol[i+1+(Nx+2)*j+1]+sol[i-1+(Nx+2)*j + 1])/(dx*dx) + (sol[i+(Nx+2)*(j+1) + 1]+sol[i+(Nx+2)*(j-1) + 1])/(dy*dy) - f[i + (Nx+2)*j + 1]);
        }
    }

    // Swap pointers
    sol.swap(solNew);
    
    //Incrementation
    //res = 
    ++l;
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
    
    for (int i = 0; i < (Nx+2)*(Ny+2); i++){
      err[i] = sol[i] - u_solution[i];
    }

    //Print
    cout << "Norm of the solution:   " << l2_norm(sol) << endl;
    cout << "Norm of expected solution:   " << l2_norm(u_solution) << endl;
    cout << "Absolute error: " << l2_norm(err) << endl;

    //File
    ofstream file;
    file.open("Jacobi.dat");
    for (int i=0; i<=Nx+2; i++){
          for (int j=0; j<=Ny+2; j++){
            file << pos_x[i] << ";" << pos_y[j] << ";" << sol[i + Nx*j] <<endl;
      }
    }
    file.close();

  /* Solution display - end */

  return 0;
}
