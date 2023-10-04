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

// Problem parameters
  double a = 1.;
  double b = 1.;
  double alpha = 0.5;
  double U_0 = 1;

  double L = 1e5;
  int Nx = 18;
  int Ny = 18;

  double dx = a/(Nx+1);
  double dy = b/(Ny+1);

  double coeff = (pow(dx,2) * pow(dy,2))/(2*(pow(dx,2) + pow(dy,2)));

// Memory allocation + Initial solution + Boundary conditions
  vector<double> f((Nx+2)*(Ny+2));
  vector<double> sol((Nx+2)*(Ny+2));
  vector<double> solNew((Nx+2)*(Ny+2));

  vector<double> u_solution((Nx+2)*(Ny+2)); //Validation

  //IC
  for (int i=0; i < (Nx+2)*(Ny+2); i++){
    sol[i] = U_0;
  }

  /*
  //LC, of problem
  for (int j = 0; j < Ny+1; j++){
    sol[0 + Nx*j] = U_0*(1+alpha*V(Nx*j,b));
  }
  */

  //f, for validation
  for (int i = 0; i<Nx+2;i++){
    for (int j=0; j<Ny+2;j++){
      f[i + (Nx+2)*j] = -2*pow(M_PI,2)*cos(M_PI*i)*sin(M_PI*j);
      u_solution[i + (Nx+2)*j] = cos(M_PI*i)*sin(M_PI*j);
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

  // Time loop
  for (int l=1; l<=L; l++){
    
    // Spatial loop
    for (int i=1; i<=Nx; i++){
        for (int j=1; j<=Ny; j++){
            solNew[i+(Nx+2)*j + 1] = coeff * ((sol[i+1+(Nx+2)*j +1]+sol[i-1+(Nx+2)*j + 1])/pow(dx,2) + (sol[i+(Nx+2)*(j+1) + 1]+sol[i+(Nx+2)*(j-1) + 1])/pow(dy,2) - f[i + (Nx+2)*j]);
        }
    }

    // Swap pointers
    sol.swap(solNew);
  }


  // Print solution

  //Def of pos
  vector<double> pos_x((Nx+2));
  vector<double> pos_y((Ny+2));

  for (int i = 0; i<Nx+2; ++i){
    pos_x[i] = i*dx;
  }

  for (int j=0; j<Ny+2; ++j){
    pos_y[j] = j*dy;
  }

  //Print
  ofstream file;
  file.open("Jacobi.txt");
  for (int i=0; i<=Nx+2; i++){
        for (int j=0; j<=Ny+2; j++){
          file << pos_x[i] << ";" << pos_y[j] << ";" << sol[i + Nx*j] <<endl;
    }
  }
  file.close();


  // Show results
  cout << "Final result:   " << l2_norm(sol) << endl;
  cout << "Norm of sol_th:   " << l2_norm(u_solution) << endl;
    //Absolute error 
    vector<double> err((Nx+2)*(Ny+2));
    for (int i = 0; i < Nx*Ny; i++){
      err[i] = sol[i] - u_solution[i];
    }
  cout << "Absolute error: " << l2_norm(err) << endl;

  return 0;
}