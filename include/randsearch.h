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
#ifndef _RANDSEARCH_H_
#define _RANDSEARCH_H_

using namespace std;

namespace nptm
{
  class RandSearch
  {

    public:
      explicit RandSearch(int num_dim, double* ptr_constrains, int type_sol, int max_iter, double sigma, int num_workers);
      void setNumDim(int num_dim);
      int getNumDim();

      void setProblemType(int type_sol);
      int getProblemType();

      void setMaxIter(int max_iter);
      int getMaxIter();

      void setConstrains(double *ptr_constrains);
      void displayConstrains();

      void setSigma(double _sigma);
      double getSigma();

      void setNumWorkers(int _num_workers);
      int getNumWorkers();

      double computeCost(double* x);
      void initSolution();
      bool checkRangeSolution(double* x);
      void updateSolution();
      void checkBestWorker(int idx_worker, double *x1, double c1, double* x2, double c2);
      void displaySolution(int iter);
      //void retrieveBestSolution();
      void optimizeSolution();

      std::function<double(double*)> cost_func;
      vector<double> constrains;
      vector<double> dim_ranges;
      vector<vector<double>> solution;
      vector<double> def_solution;
      vector<double> fitness;


    private:
      int type_sol;
      int num_dim;
      int max_iter;
      double sigma;
      int num_workers;
  };
}

#endif
