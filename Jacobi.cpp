// Compilation:
//   g++ Jacobi.cpp
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
        accum += u[i]*u[i];
    }
    return sqrt(accum/u.size());
}

double V(double y, double b){
    return 1 - cos(2*M_PI*y/b);
}

int main(int argc, char* argv[]){

//Algorithm parameters
  double L = 1e5;
  double epsilon = 1e-5;

// Problem parameters
  double a = 1.;
  double b = 1.;
  double alpha = 0.5;
  double U_0 = 0;
  int Nx = 30;
  int Ny = 30;

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
  

  /* Validation process - begin*/

  
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

  /* Validation process - end */


  /* Time loop - begin*/

  //Initiallization
  int l = 0;
  vector<double> res((Nx+2)*(Ny+2));
  double n_res = MAXFLOAT;
  double n_res_0 = MAXFLOAT;

  // Loop

    // l = 0
  while (l < 1){
    // Spatial loop
    for (int i=1; i<=Nx; i++){
        for (int j=1; j<=Ny; j++){
            solNew[i+(Nx+2)*j + 1] = coeff * ((sol[i+1+(Nx+2)*j+1]+sol[i-1+(Nx+2)*j + 1])/(dx*dx) + (sol[i+(Nx+2)*(j+1) + 1]+sol[i+(Nx+2)*(j-1) + 1])/(dy*dy) - f[i + (Nx+2)*j + 1]);     
        }
    }
    // Swap pointers
    sol.swap(solNew);
    
    //Residu
    for (int i = 0; i < (Nx+2)*(Ny+2); i++){
      res[i] = sol[i] - solNew[i];
    }
    n_res_0 = l2_norm(res);

    //Incrementation
    l++;
  }

    // l > 1
  while((l<=L)&&(n_res/n_res_0 > epsilon)){
    
    // Spatial loop
    for (int i=1; i<=Nx; i++){
        for (int j=1; j<=Ny; j++){
            solNew[i+(Nx+2)*j + 1] = coeff * ((sol[i+1+(Nx+2)*j+1]+sol[i-1+(Nx+2)*j + 1])/(dx*dx) + (sol[i+(Nx+2)*(j+1) + 1]+sol[i+(Nx+2)*(j-1) + 1])/(dy*dy) - f[i + (Nx+2)*j + 1]);     
        }
    }

    // Swap pointers
    sol.swap(solNew);
    
    //Residu
    for (int i = 0; i < (Nx+2)*(Ny+2); i++){
      res[i] = sol[i] - solNew[i];
    }
    n_res = l2_norm(res);
    
    //Incrementation
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
      err[i] = sol[i] - sol_th[i];
    }
    

    //Print
    cout << "***** Printing the solution *****" <<endl;

    //finding out why the program stopped
    if (l > L){
      cout << "Program stopped beacause l = " << l << "> L" << endl;
    }

    else {
      cout << "Program stopped because the residual value = " << n_res/n_res_0 << "< "<< epsilon << endl;
    }

    cout << "- Absolute error: " << l2_norm(err) << endl;

    cout << "- Norm of the residual value: " << n_res/n_res_0 << endl;

    cout << "********************************" <<endl;

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
