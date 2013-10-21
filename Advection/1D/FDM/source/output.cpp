#include <output.h>
#include <norms.h>

#include <Eigen/Dense>

#include <iostream>
#include <iomanip>
#include <ostream>
#include <fstream>

using namespace Eigen;
using namespace std;

void displayField (VectorXd u) {
  cout << u << endl;
}

void displayField (VectorXd u,
                   ostream & stream) {
  stream << u << endl;
}

void outputStats (const int N,
                  const double h,
                  const double dt,
                  const int stride,
                  const double T) {
  std::cerr << "==> Run statistics:" << endl
            << "====> Subdivisions (N)       = " << N << endl
            << "====> Timestep (dt)          = " << dt << endl
            << "====> Output stride (stride) = " << stride << endl
            << "====> End time (T)           = " << T << endl << endl;
}

ofstream & openTeXDoc (string filename) {
  ofstream * stream_ptr = new ofstream;
  ofstream & stream     = *stream_ptr;
  stream.open (filename.c_str(), ofstream::out | ofstream::trunc);

  stream << "\\documentclass[12pt]{article}" << endl
         << "\\begin{document}" << endl
         << "\t\\begin{tabular}{l|c|c|c|c|c|c}" << endl
         << "\t\tN&Linf&Convergence&L1&Convergence&L2&Convergence\\\\" << endl
         << "\t\t&&Rate&&Rate&&Rate\\\\" << endl;

  return stream;
}

void makeTeXRow (ofstream & stream,
                 VectorXd error) {
  static double om;
  static double oo;
  static double ot;

  stream << "\t\t\\hline" << endl;

  double nm = maxNorm (error);
  double no = oneNorm (error);
  double nt = twoNorm (error);

  int N = error.rows();

  stream << "\t\t" << N
         << "&" << scientific << setprecision (6) << nm << fixed
         << "&" << setprecision (2) << ((om != 0) ? log (om / nm) / log (2.0) : 0)
         << "&" << scientific << setprecision (6) << no << fixed
         << "&" << setprecision (2) << ((oo != 0) ? log (oo / no) / log (2.0) : 0)
         << "&" << scientific << setprecision (6) << nt << fixed
         << "&" << setprecision (2) << ((ot != 0) ? log (ot / nt) / log (2.0) : 0)
         << "\\\\" << endl;
  om = nm; oo = no; ot = nt;
}

void closeTeXDoc (std::ofstream & stream) {
  stream << "\t\\end{tabular}" << endl
         << "\\end{document}" << endl;
  stream.close();
}
