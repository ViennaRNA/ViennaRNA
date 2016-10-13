#include <stdio.h>
#include <iostream>
#include <g2_PS.h>
#include <g2.h>

/* A very simple example to demonstrate g2 in a C++ environment */


class Circle
{
public:
  Circle(int d, double x, double y, double r)
    : _d(d), _x(x), _y(y), _r(r)
  {
    g2_circle(d, x, y, r);
    g2_flush(d);
  }

  void Fill()
  {
    g2_filled_circle(_d, _x, _y, _r);
    g2_flush(_d);
  }
private:
  int _d;
  double _x, _y, _r;
};


int main(int argc, char *argv[])
{
    int d;
    d=g2_open_PS("demo_cpp.ps", g2_A4, g2_PS_port);
    Circle c(d, 150, 150, 25);
    c.Fill();
    Circle c2(d, 150, 175, 25);
    g2_close(d);
    std::cout << "demo_cpp.ps generated\n";
    return 0;
}
