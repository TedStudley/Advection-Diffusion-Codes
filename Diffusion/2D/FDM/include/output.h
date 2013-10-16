#pragma once

#include <Eigen/Dense>

#include <fstream>


void displayField (Eigen::VectorXd);

void displayField (Eigen::VectorXd,
                   std::ofstream &);

std::ofstream & openTeXDoc (std::string);

void makeTeXRow (std::ofstream &,
                 Eigen::VectorXd);

void closeTeXDoc (std::ofstream &);
