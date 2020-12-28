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

#include "simulation.h"

class Window : public QWidget
{
    Q_OBJECT

 public:
    Window(gridDim dim, float3 B, float J, float beta,
	   simulationType simulationType, bool heat, QWidget *parent = 0);

    void paintEvent(QPaintEvent *);

 private:
    Models models;
    std::unique_ptr<Grid<short int>> grid;
};

#endif
