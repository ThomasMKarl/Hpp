#ifndef WINDOW_H
#define WINDOW_H

#include <QWidget>
#include <QTimer>
#include <string>
#include <QLineEdit>
#include <QLabel>
#include <QGridLayout>
#include <iostream>
#include <QPainter>

#include "qt_application/app.h"
#include "simulation.h"

enum SimulationType {metropolis, heatbath, wolff};

typedef struct {
  size_t x;
  size_t y;
}gridSizes;

class Window : public QWidget
{
    Q_OBJECT

 public:
    Window(gridSizes sizes,
	   float3 B, float J, float beta,
	   SimulationType simulationType, bool heat, QWidget *parent = 0);

    void paintEvent(QPaintEvent *);

 private:
    HB::Models models;
    std::unique_ptr<HB::Grid<short int>> grid;
};

#endif
