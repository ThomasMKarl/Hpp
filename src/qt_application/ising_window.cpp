#include "qt_application/ising_window.h"

Window::Window(gridSizes sizes, float3 B, float J, float beta,
	       SimulationType simulationType, bool heat, QWidget *parent) : QWidget(parent)
{
    dim3 gridSizes3d{sizes.x,sizes.y,1};
    short int spatialDimension = 2;
    grid = std::make_unique<HB::Grid<short int>>(spatialDimension, gridSizes3d);
    grid.calcNeighbourTable();
    if(heat) grid.hotStart();
    else     grid.coldStart();

    if(simulationType == metropolis) HB::MetropolisSimulationQt simStrategy{grid,beta};
    //if(simulationType == wolff)      simStrategy(grid,beta,steps);
    //if(simulationType == heatbath)   simStrategy(grid,beta,steps);

    models.push_back(std::make_unique<HB::Model>>(J, B, simStrategy));
}

void Window::paintEvent(QPaintEvent *) //override paintEvent function
{
    if(simulationType == metropolis) HB::simulate(models);

    /*if(wmc != nullptr) //if wolff algorithm
    {
        wmc->MC_sweep();

        QPen WhitePen(Qt::white, PenSize, Qt::SolidLine),
             BlackPen(Qt::black, PenSize, Qt::SolidLine),
               RedPen(Qt::red,   PenSize, Qt::SolidLine);

        for (unsigned i = 0; i < wmc->XSize(); i++)
        {

              for (unsigned j = 0; j < wmc->YSize(); j++)
              {

                  if(wmc->compute_spin(wmc->XSize()*j+i) ==  1) qp.setPen(WhitePen);
                  if(wmc->compute_spin(wmc->XSize()*j+i) == -1)	qp.setPen(BlackPen);
                  if(wmc->compute_spin(wmc->XSize()*j+i) == -2 || wmc->compute_spin(wmc->XSize()*j+i) == 2) //2 and -2 refer to spins of a cluster
                      qp.setPen(RedPen);

                  qp.drawPoint(offset * i, offset * j); //draw white point if spin 1, black if -1, draw points of cluster (2 or -2) red

              }
        }
	}*/
}
