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

#include <iostream>
#include <stdio.h>
#include <cmath>
#include <random>
#include <chrono>
#include <string>
#include <bits/stdc++.h>

#include "randsearch.h"

using namespace std;

namespace nptm
{
  RandSearch::RandSearch(int _num_dim, double* ptr_constrains, int _type_sol, int _max_iter, double _sigma, int _num_workers)
  {
    setNumDim(_num_dim);
    setConstrains(ptr_constrains);
    setProblemType(_type_sol);
    setMaxIter(_max_iter);
    setSigma(_sigma);
    setNumWorkers(_num_workers);
    cout << "random search problem created! " << endl;
  }

  void RandSearch::setNumDim(int _num_dim) { num_dim = _num_dim; }
  int RandSearch::getNumDim() { return num_dim; }

  void RandSearch::setProblemType(int _type_sol) { type_sol = _type_sol; }
  int RandSearch::getProblemType() { return type_sol; }

  void RandSearch::setMaxIter(int _max_iter) { max_iter = _max_iter; }
  int RandSearch::getMaxIter() { return max_iter; }

  void RandSearch::setSigma(double _sigma) { sigma = _sigma; }
  double RandSearch::getSigma() { return sigma; }

  void RandSearch::setNumWorkers(int _num_workers) { num_workers = _num_workers; }
  int RandSearch::getNumWorkers() { return num_workers; }

  void RandSearch::setConstrains(double* ptr_constrains)
  {
    for(unsigned int i=0; i < getNumDim()*2; i++)
    {
      constrains.push_back(ptr_constrains[i]);
      cout << "constrains[" << i << "] = " << constrains.at(i) << "\t" << endl;
    }

    //cout << "constrains size = " << constrains.size() << endl;

    for(unsigned int i=0; i < getNumDim(); i++)
    {
      //double range_current_dim = ;
      dim_ranges.push_back(constrains.at(2*i+1) - constrains.at(2*i));
      cout << "range in dimension " << i << ": " << dim_ranges.at(i) << "\t" << endl;
    }
  }

  void RandSearch::displayConstrains()
  {
    for(int i=0; i < getNumDim(); i++)
    {
      cout << "dimension " << i << ":" << endl;
      cout << "x" << i << "_min = " << constrains.at(2*i) << "\t";
      cout << "x" << i << "_max = " << constrains.at(2*i+1) << endl;
    }
  }


  double RandSearch::computeCost(double* x)
  {
    return cost_func(x);
  }

  void RandSearch::initSolution()
  {
    for(unsigned int i = 0; i < getNumWorkers(); i++)
    {
      // init definitive solution
      def_solution.push_back(0.0);

      // init fitness for each worker
      fitness.push_back(-100000.0);

      // init the solution for each worker
      vector<double> temp;
      for(unsigned int j = 0; j < getNumDim(); j++)
      {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        uniform_real_distribution<double> distribution(0.0, 1.0);
        temp.push_back(distribution(generator)*dim_ranges.at(j) + constrains.at(2*j));
      }
      solution.push_back(temp);
    }
  }

  bool RandSearch::checkRangeSolution(double* x)
  {
    for(unsigned int i=0; i < getNumDim(); i++)
    {
      // test if solution in dimension i is smaller than the minimum
      // or if the solution in dimension i is larger than the maximum
      if((x[i] < constrains.at(2*i)) || (x[i] > constrains.at(2*i+1)))
      {
        return false;
      }
    }
    return true;
  }

  void RandSearch::updateSolution()
  {
    for(unsigned int i = 0; i < getNumWorkers(); i++)
    {
      for(unsigned int j=0; j < getNumDim(); j++)
      {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        uniform_real_distribution<double> distribution(0.0, 1.0);
        solution[i][j] = solution[i][j] + distribution(generator)*getSigma();
      }
    }
  }


  void RandSearch::checkBestWorker(int idx_worker, double* x_1, double c1, double* x_2, double c2)
  {
    if(c2 > c1)
    {
      fitness[idx_worker] = c1;
      for(unsigned int i=0; i < getNumDim(); i++)
      {
        solution[idx_worker][i] = x_1[i];
        //cout << "x[" << i << "] = " << x[i] << endl;
      }
    } else {
      fitness[idx_worker] = c2;
      for(unsigned int i=0; i < getNumDim(); i++)
      {
        solution[idx_worker][i] = x_2[i];
        //cout << "x[" << i << "] = " << x[i] << endl;
      }
    }
  }

  void RandSearch::displaySolution(int _iter)
  {
    // first display the current iteration plus the columns of the data
    cout << "Iteration " << _iter << "  #################################" << endl;

    cout << "Worker Id" << "\t" << "Fitness" << "\t";

    for(unsigned int j=0; j < getNumDim(); j++)
    {
      cout << "x" << j << "\t";
    }
    cout << endl;


    // now display the current values for the iteration
    for(unsigned int i=0; i < getNumWorkers(); i++)
    {
      cout << i << "\t" << fitness[i] << "\t";

      for(unsigned int j=0; j < getNumDim(); j++)
      {
        cout << solution[i][j] << "\t";
      }
      cout << endl;
    }
  }



  /*
  void RandSearch::retrieveBestSolution()
  {
      for(unsigned int i=0; i < getNumWorkers(); i++)
      {
        for
      }
  }
  */

  void RandSearch::optimizeSolution()
  {
    // INITIALIZATION /////////////////////////////////////////////////////////
    //cout << "initial solution = " << endl;
    initSolution();
    double x_1[getNumWorkers()][getNumDim()] = {0.0};
    double x_2[getNumWorkers()][getNumDim()] = {0.0};
    double cost_1[getNumWorkers()] ={0.0};
    double cost_2[getNumWorkers()] = {0.0};
    for(unsigned int i=0; i < getNumWorkers(); i++)
    {
      for(unsigned int j=0; j < getNumDim(); j++)
      {
        x_1[i][j] = solution[i][j];
        //cout << "x_1[" << i << "] = " << x_1[i] << endl;
      }
    }


    // SEARCH LOOP ////////////////////////////////////////////////////////////
    for(unsigned int i=0; i < max_iter; i++)
    {
      for(unsigned int j=0; j < getNumWorkers(); j++)
      {
        double current_worker_solution[getNumDim()];
        for(unsigned int k=0; k < getNumDim(); k++)
        {
          current_worker_solution[k] = x_1[j][k];
        }
        if(checkRangeSolution(current_worker_solution))
        {
          cost_1[j] = computeCost(current_worker_solution);
        } else {
          cost_1[j] = 10000;
        }
      }


      // update solution class member
      for(unsigned int j=0; j < getNumWorkers(); j++)
      {
        for(unsigned k=0; k < getNumDim(); k++)
        {
          solution[j][k] = x_1[j][k];
          //cout << "x[" << i << "] = " << x[i] << endl;
        }
      }

      updateSolution();

      for(unsigned int j=0; j < getNumWorkers(); j++)
      {
        for(unsigned k=0; k < getNumDim(); k++)
        {
          x_2[j][k] = solution[j][k];
          //cout << "x_2[" << j << "] = " << x_2[j] << endl;
        }
      }


      for(unsigned int j=0; j < getNumWorkers(); j++)
      {
        double current_worker_solution[getNumDim()];
        for(unsigned int k=0; k < getNumDim(); k++)
        {
          current_worker_solution[k] = x_2[j][k];
        }
        if(checkRangeSolution(current_worker_solution))
        {
          cost_2[j] = computeCost(current_worker_solution);
        } else {
          cost_2[j] = 10000;
        }
      }

      // UPDATE SOLUTION //////////////////////////////////////////////////////
      double current_x_1[getNumDim()] = {0.0};
      double current_x_2[getNumDim()] = {0.0};
      for(unsigned int j=0; j < getNumWorkers(); j++)
      {
        for(unsigned k=0; k < getNumDim(); k++)
        {
          current_x_1[k] = x_1[j][k];
          current_x_2[k] = x_2[j][k];
        }
        checkBestWorker(j, current_x_1, cost_1[j], current_x_2, cost_2[j]);
      }

      // DISPLAY THE SOLUTIONS FOR EACH ITERATION
      displaySolution(i);

      // RETRIEVE BEST SOLUTION PER ITERATION

    } // END OF SEARCH LOOP
  }
}
