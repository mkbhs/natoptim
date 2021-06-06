#include <iostream>
#include <cmath>
#include <random>
#include <ctime>

#include <bits/stdc++.h>

#include "randsearch.h"

using namespace std;

namespace nptm
{
  RandSearch::RandSearch()
  {
    x = 0;
    y = 0;
  }

  int RandSearch::sum(const int _x, const int _y)
  {
    x = _x;
    y = _y;
    return x+y;
  }

  /*
  int RandSearch::computeCostFunction(int x, int y, std::function<int(int, int)> func)
  {
    return func(x,y);
  }
  */
  int RandSearch::computeCostFunction()
  {
    return func(x,y);
  }
  
}
