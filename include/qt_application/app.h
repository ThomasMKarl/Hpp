#ifndef APP_H
#define APP_H

#include <QWidget>
#include <QPushButton>
#include <QLineEdit>
#include <QGridLayout>
#include <QLabel>
#include <QTimer>
#include <string>
#include <QApplication>
#include <QPainter>
#include <iostream>

#include "qt_application/ising_window.h"
  
class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent = 0);

private slots:
    void metroButtonClicked();
    void    hbButtonClicked();
    void wolffButtonClicked();
    void   hotButtonClicked();
    void  coldButtonClicked();

    void start_ising_window(SimulationType simulationType);

    template<typename T>
    bool isCharacter(T value, QLineEdit edit, QString errorMessage);
    template<typename T>
    bool isEmpty(T value, QLineEdit &edit, QString errorMessage);
    void connectTimerWithWindow(QTimer &timer, Window &window);

private:
    QLineEdit *displayT, *displayJ, *displayBx, *displayBy, *displayBz;
    QLabel *T, *kbT, *J, *X, *Y, *B, *displayX, *displayY;
    QPushButton *hbButton, *metroButton, *wolffButton, *hotButton, *coldButton;
    QGridLayout mLayout;

    bool heat = true;
};

#endif // APP_H
