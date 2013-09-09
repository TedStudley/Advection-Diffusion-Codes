#include <initial-conditions.h>
#include <diffusion.h>
#include <output.h>
#include <norms.h>

#include <Eigen/Dense>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>

using namespace Eigen;
using namespace std;


int main() {
  const double mu      = 0.4;
  const double kappa   = 1.0;
  const double t0      = 0.0001;
  const int k = 1;
  double t;

  ofstream outstream;
  string output_dir = "output/";
  string file_prefix = "FE-";
  stringstream filename;

  for (int N = 16; N <= 2048; N*=2) {
    VectorXd u (N);
    VectorXd u1 (N);
    VectorXd utemp (N);

    double h  = 1.0 / (N + 1);
    double delta_t = mu * h / kappa;

    squareWave (u1);
    fourierSquare (u, kappa, t0);

    t = t0;

    for (int i = 0; t < 0.05; ++i) {
      crankNicolson (u, delta_t, h, kappa);
      t += delta_t;

      if (i % 5 == 0) {
        filename << output_dir << file_prefix 
                 << setw (6) << setfill('0') << N << "-" 
                 << setw (6) << setfill ('0') << i << ".dat";
        outstream.open(filename.str().c_str());

        displayField(u, outstream);
        sineWave(utemp, k, kappa, t);
        displayField(utemp, outstream);

        filename.str(std::string());
        outstream.close();
      }

    }

    fourierSquare (u1, kappa, t);


    VectorXd error = (u - u1);

    cout << N << " " << maxNorm (u - u1) 
         << " " << sqrt(error.array().square().sum()) / N << endl;

  }

  return 0;
}
