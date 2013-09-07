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
  const int k = 1;
  double t;

  ofstream outstream;
  string output_dir = "output/";
  string file_prefix = "FE-";
  stringstream filename;

  for (int N = 16; N <= 512; N*=2) {
    VectorXd u (N);
    VectorXd u1 (N);

    double h  = 1.0 / (N + 1);
    double delta_t = mu * h / kappa;

    sineWave (u, N, k);

    t = 0;

    for (int i = 0; t < 0.05; ++i) {
      crankNicolson (u, N, mu * N);
      t += delta_t;

      if (i % 5 == 0) {
        filename << output_dir << file_prefix 
                 << setw (6) << setfill('0') << N << "-" 
                 << setw (6) << setfill ('0') << i << ".dat";
        outstream.open(filename.str().c_str());

        displayField(u, N, outstream);
        sineWave(u1, N, k, kappa, t);
        displayField(u1, N, outstream);

        filename.str(std::string());
        outstream.close();
      }

    }

    sineWave (u1, N, k, kappa, t);


    VectorXd error = (u - u1);

    cout << N << " " << maxNorm (u - u1) 
      << " " << sqrt(error.array().square().sum()) / N << endl;

  }

  return 0;
}
