
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
