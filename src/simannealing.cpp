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

#include "simannealing.h"

using namespace std;

namespace nptm
{
  SimAnnealing::SimAnnealing(int _num_dim, double* ptr_constrains, int _max_iter, double _i_temp, double _f_temp, double _beta, double _sigma, int _num_workers)
  {
    setNumDim(_num_dim);
    setConstrains(ptr_constrains);
    setMaxIter(_max_iter);
    setInitialTemp(_i_temp);
    setFinalTemp(_f_temp);
    setBeta(_beta);
    setSigma(_sigma);
    setNumWorkers(_num_workers);
  }

  void SimAnnealing::setNumDim(int _num_dim) { num_dim = _num_dim; }
  int SimAnnealing::getNumDim() { return num_dim; }

  void SimAnnealing::setMaxIter(int _max_iter) { max_iter = _max_iter; }
  int SimAnnealing::getMaxIter() { return max_iter; }

  void SimAnnealing::setInitialTemp(double _i_temp) { i_temp = _i_temp; }
  double SimAnnealing::getInitialTemp() { return i_temp; }

  void SimAnnealing::setFinalTemp(double _f_temp) { f_temp = _f_temp; }
  double SimAnnealing::getFinalTemp() { return f_temp; }

  void SimAnnealing::setCurrentTemp(double _c_temp) { c_temp = _c_temp; }
  double SimAnnealing::getCurrentTemp() { return c_temp; }

  void SimAnnealing::setSigma(double _sigma) { sigma = _sigma; }
  double SimAnnealing::getSigma() { return sigma; }

  void SimAnnealing::setBeta(double _beta) { beta = _beta; }
  double SimAnnealing::getBeta() { return beta; }

  void SimAnnealing::setNumWorkers(int _num_workers) { num_workers = _num_workers; }
  int SimAnnealing::getNumWorkers() { return num_workers; }

  void SimAnnealing::setDefFitness(double _def_fitness) { def_fitness = _def_fitness; }
  double SimAnnealing::getDefFitness() { return def_fitness; }



  void SimAnnealing::setConstrains(double* ptr_constrains)
  {
    for(unsigned int i=0; i < getNumDim()*2; i++)
    {
      constrains.push_back(ptr_constrains[i]);
      cout << "constrains[" << i << "] = " << constrains.at(i) << "\t" << endl;
    }

    for(unsigned int i=0; i < getNumDim(); i++)
    {
      dim_ranges.push_back(constrains.at(2*i+1) - constrains.at(2*i));
      cout << "range in dimension " << i << ": " << dim_ranges.at(i) << "\t" << endl;
    }
  }


  void SimAnnealing::displayConstrains()
  {
    for(int i=0; i < getNumDim(); i++)
    {
      cout << "dimension " << i << ":" << endl;
      cout << "x" << i << "_min = " << constrains.at(2*i) << "\t";
      cout << "x" << i << "_max = " << constrains.at(2*i+1) << endl;
    }
  }


  double SimAnnealing::computeCost(double* x)
  {
    return cost_func(x);
  }



  void SimAnnealing::initSolution()
  {
    for(unsigned int i = 0; i < getNumWorkers(); i++)
    {
      // init definitive solution
      def_solution.push_back(0.0);

      // init fitness for each worker
      fitness.push_back(100000.0);

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


  void SimAnnealing::updateSolution()
  {
    for(unsigned int i = 0; i < getNumWorkers(); i++)
    {
      for(unsigned int j=0; j < getNumDim(); j++)
      {
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        uniform_real_distribution<double> distribution(0.0, 1.0);
        solution[i][j] = solution[i][j] + distribution(generator)*getSigma(); //*getCurrentTemp();
      }
    }
  }


  bool SimAnnealing::checkRangeSolution(double* x)
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

  void SimAnnealing::displaySolution(int _iter)
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


  void SimAnnealing::retrieveBestSolutionGlobal()
  {
    int idx_min = 0;
    double min_cost = 10000000.0;
    double current_cost = 0.0;
    for(unsigned int i=0; i < getNumWorkers(); i++)
    {
      double best_solution[getNumDim()] = {0.0};
      for(unsigned int j=0; j < getNumDim(); j++)
      {
        best_solution[j] = solution[i][j];
      }
      // check the cost
      current_cost = computeCost(best_solution);
      if(current_cost < min_cost)
      {
        min_cost = current_cost;
        idx_min = i;
      }
    }
    // update the definitive solution
    for(unsigned int i=0; i < getNumDim(); i++)
    {
      def_solution[i] = solution[idx_min][i];
    }
    setDefFitness(min_cost);
  }


  void SimAnnealing::displayGlobalSolution()
  {
    cout << "Best Solution" << endl;
    for(unsigned int i=0; i < getNumDim(); i++)
    {
      cout << "x" << i << "\t";
    }
    cout << endl;
    for(unsigned int i=0; i < getNumDim(); i++)
    {
      cout << def_solution[i] << "\t";
    }
    cout << endl;
    cout << "with a cost = " << getDefFitness() << endl;
  }


  void SimAnnealing::optimizeSolution()
  {
    // INITIALIZATION /////////////////////////////////////////////////////////
    initSolution();
    //initBestSolution();
    setCurrentTemp(getInitialTemp());

    double x_1[getNumWorkers()][getNumDim()] = {0.0};
    double x_2[getNumWorkers()][getNumDim()] = {0.0};

    for(unsigned int i=0; i < getNumWorkers(); i++)
    {
      for(unsigned int j=0; j < getNumDim(); j++)
      {
        x_1[i][j] = solution[i][j];
      }
    }

    // OPTIMIZATION LOOP
    for(unsigned int i=0; i < getMaxIter(); i++)
    {
      double E_old[getNumWorkers()] = {0.0};
      for(unsigned int j=0; j < getNumWorkers(); j++)
      {
        double current_worker_solution[getNumDim()];
        for(unsigned int k=0; k < getNumDim(); k++)
        {
          current_worker_solution[k] = x_1[j][k];
        }
        E_old[j] = computeCost(current_worker_solution);
      }


      updateSolution();
      for(unsigned int j=0; j < getNumWorkers(); j++)
      {
        for(unsigned k=0; k < getNumDim(); k++)
        {
          x_2[j][k] = solution[j][k];
        }
      }

      double E_new[getNumWorkers()] = {0.0};
      for(unsigned int j=0; j < getNumWorkers(); j++)
      {
        double current_worker_solution[getNumDim()];
        for(unsigned int k=0; k < getNumDim(); k++)
        {
          current_worker_solution[k] = x_2[j][k];
        }
        E_new[j] = computeCost(current_worker_solution);
      }

      // compute deltas
      double delta[getNumWorkers()] = {0.0};
      for(unsigned int j=0; j < getNumWorkers(); j++)
      {
        // if cost is smaller than second solution, accept it directly
        delta[j] = E_new[j] - E_old[j];
        if(delta[j] < 0)
        {
          for(unsigned int k=0; k < getNumDim(); k++)
          {
            x_1[j][k] = x_2[j][k];
            solution[j][k] = x_1[j][k];
          }
          fitness[j] = E_new[j];
        } else {
          for(unsigned int k=0; k < getNumDim(); k++)
          {
            solution[j][k] = x_1[j][k];
          }
          fitness[j] = E_old[j];
        }

        // if the cost is not smaller, check the acceptance probability
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        default_random_engine generator(seed);
        uniform_real_distribution<double> distribution(0.0, 1.0);
        if((delta[j] >= 0) && (exp(delta[j]/getCurrentTemp()) < distribution(generator)))
        {
          for(unsigned int k=0; k < getNumDim(); k++)
          {
            x_1[j][k] = x_2[j][k];
            solution[j][k] = x_1[j][k];
          }
          fitness[j] = E_new[j];
        } else {
          for(unsigned int k=0; k < getNumDim(); k++)
          {
            solution[j][k] = x_1[j][k];
          }
          fitness[j] = E_old[j];
        }
      }

      // Display the solutions for each iteration
      //displaySolution(i);


      // UPDATE THE CURRENT TEMPERATURE
      double current_temp = getBeta()*getCurrentTemp();

      if(current_temp <= getFinalTemp())
      {
        setCurrentTemp(getFinalTemp());
      } else {
        setCurrentTemp(current_temp);
      }
      //cout << "current teperature = " << getCurrentTemp() << endl;

    } // END OF OPTIMIZATION LOOP

    // Retrieve the best solution among all workers
    retrieveBestSolutionGlobal();

    // display best solution
    displayGlobalSolution();

  }

}
