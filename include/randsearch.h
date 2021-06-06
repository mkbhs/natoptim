#pragma once
#ifndef _RANDSEARCH_H_
#define _RANDSEARCH_H_

using namespace std;

namespace nptm
{
  class RandSearch
  {

    public:
      explicit RandSearch();
      int sum(int x, int y);
      //int computeCostFunction(int x, int y,std::function<int(int, int)> func);
      int computeCostFunction();

      int x, y;
      std::function<int(int, int)> func;

    //private:
    //  int x, y;

  };
}

#endif
