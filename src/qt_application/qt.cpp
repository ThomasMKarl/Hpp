/*--------------------------------------------

PROJECT: Visualisation of Ising Model with Qt

Thomas Karl

--------------------------------------------*/

#include "qt_application/app.h"
#include <QApplication>

int main(int argc, char *argv[])
{
    QApplication a(argc, argv);
    Widget w;
    w.show();

    return a.exec();
}
