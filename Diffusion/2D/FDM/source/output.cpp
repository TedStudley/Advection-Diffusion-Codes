#include <output.h>
#include <norms.h>

#include <Eigen/Dense>

#include <iostream>
#include <iomanip>
#include <fstream>


void displayField (Eigen::VectorXd u) {
  const int N = sqrt (u.rows ());

  for (int i = 0; i < N; ++i)
    std::cout << u.segment (i * N, N).transpose () << std::endl; 
}

void displayField (Eigen::VectorXd u,
                   std::ofstream & stream) {
  const int N = sqrt (u.rows ());

  for (int i = 0; i < N; ++i)
    stream << u.segment (i * N, N).transpose () << std::endl;
}

std::ofstream & openTeXDoc (std::string filename) {
  std::ofstream * stream_ptr = new std::ofstream;
  std::ofstream & stream = *stream_ptr;
  stream.open (filename.c_str(), std::ofstream::out | std::ofstream::trunc);

  stream << "\\documentclass[12pt]{article}" << std::endl
         << "\\begin{document}" << std::endl;
  stream << "\t\\begin{tabular}{l|c|c|c|c|c|c}" << std::endl
         << "\t\tN&Max Norm&Convergence&One Norm&Convergence&Two Norm&Convergence\\\\" << std::endl
         << "\t\t&&Rate&&Rate&&Rate\\\\" << std::endl;

  return stream;
}

void makeTeXRow (std::ofstream & stream,
                 Eigen::VectorXd error) {
  static double om;
  static double oo;
  static double ot;

  stream << "\t\t\\hline" << std::endl;

  double nm = maxNorm (error);
  double no = oneNorm (error);
  double nt = twoNorm (error);

  int N = sqrt (error.rows ());

  stream << "\t\t" << N
         << "&" << std::scientific << std::setprecision (6) << nm << std::fixed
         << "&" << std::setprecision (2) << ((om != 0) ? log (om / nm) / log (2.0) : 0)
         << "&" << std::scientific << std::setprecision (6) << no << std::fixed
         << "&" << std::setprecision (2) << ((oo != 0) ? log (oo / no) / log (2.0) : 0)
         << "&" << std::scientific << std::setprecision (6) << nt << std::fixed
         << "&" << std::setprecision (2) << ((ot != 0) ? log (ot / nt) / log (2.0) : 0)
         << "\\\\" << std::endl;
  om = nm; oo = no; ot = nt;
}

void closeTeXDoc (std::ofstream & stream) {
  stream << "\t\\end{tabular}" << std::endl
         << "\\end{document}" << std::endl;
  stream.close ();
}
