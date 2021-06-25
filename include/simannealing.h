// MIT License
//
// Copyright (c) 2021 mkbhs
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

#pragma once
#ifndef _SIMANNEALING_H_
#define _SIMANNEALING_H_

using namespace std;

namespace nptm
{
  class SimAnnealing
  {
    public:
      explicit SimAnnealing(int num_dim, double* ptr_constrains, int max_iter, double i_temp, double f_temp, double beta, double sigma, int num_workers);

      void optimizeSolution();

      std::function<double(double*)> cost_func;
      vector<double> constrains;
      vector<double> dim_ranges;
      vector<vector<double>> solution;
      vector<double> def_solution;
      vector<double> fitness;

    private:

      void setNumDim(int num_dim);
      int getNumDim();

      void setMaxIter(int max_iter);
      int getMaxIter();

      void setConstrains(double *ptr_constrains);
      void displayConstrains();

      void setInitialTemp(double i_temp);
      double getInitialTemp();

      void setFinalTemp(double f_temp);
      double getFinalTemp();

      void setCurrentTemp(double c_temp);
      double getCurrentTemp();

      void setSigma(double sigma);
      double getSigma();

      void setBeta(double beta);
      double getBeta();

      void setNumWorkers(int _num_workers);
      int getNumWorkers();

      void setDefFitness(double _def_fitness);
      double getDefFitness();

      double computeCost(double* x);
      void initSolution();
      bool checkRangeSolution(double* x);
      void updateSolution();
      void checkBestWorker(int idx_worker, double *x1, double c1, double* x2, double c2);
      void displaySolution(int iter);
      void retrieveBestSolutionGlobal();
      void displayGlobalSolution();


      int num_dim;
      int max_iter;
      double i_temp;
      double f_temp;
      double c_temp;
      double sigma;
      double beta;
      int num_workers;
      double def_fitness;
  };

}
#endif
