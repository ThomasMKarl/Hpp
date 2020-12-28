#include "qt_application/app.h"
#include "qt_application/ising_window.h"


Widget::Widget(QWidget *parent): QWidget(parent)
{
        QGridLayout *layout=new QGridLayout;

        metroButton = new QPushButton("Metropolis");
        connect(metroButton, &QPushButton::clicked, this, &Widget::metroButtonClicked);
        layout->addWidget(metroButton,4,0);

        hbButton    = new QPushButton("Heatbath");
        connect(hbButton, &QPushButton::clicked, this, &Widget::hbButtonClicked);
        layout->addWidget(hbButton,4,1);

        wolffButton = new QPushButton("Wolff");
        connect(wolffButton, &QPushButton::clicked, this, &Widget::wolffButtonClicked);
        layout->addWidget(wolffButton,4,2);

        hotButton = new QPushButton("hot start");
        connect(hotButton, &QPushButton::clicked, this, &Widget::hotButtonClicked);
        layout->addWidget(hotButton,4,3);
        hotButton->setCheckable(true);
        hotButton->setChecked(true);

        coldButton = new QPushButton("cold start");
        connect(coldButton, &QPushButton::clicked, this, &Widget::coldButtonClicked);
        layout->addWidget(coldButton,4,4);
        coldButton->setCheckable(true);

        //set default values in displays
        displayT = new QLineEdit(this);
        T   = new QLabel("Temperature:", this);
        displayT->setText("0.0001");
        kbT = new QLabel("as a multible of J/k_b", this);
	
        displayX = new QLineEdit(this);
        X   = new QLabel("grid size in x-direction:", this);
        displayY = new QLineEdit(this);
        Y   = new QLabel("grid size in y-direction:", this);

        displayJ    = new QLineEdit(this);
        J   = new QLabel("Coupling constant J:", this);
        displayJ->setText("1.0");
	
        displayBx   = new QLineEdit(this);
        Bx   = new QLabel("External magnetic field, x direction:", this);
        displayBx->setText("0.0");
	displayBy   = new QLineEdit(this);
        By   = new QLabel("External magnetic field, y direction:", this);
        displayBy->setText("0.0");
	displayBz   = new QLineEdit(this);
        Bz   = new QLabel("External magnetic field, z direction:", this);
        displayBz->setText("0.0");



        layout->addWidget(displayT,1,1);
                layout->addWidget(T,1,0);
                layout->addWidget(kbT,1,2);
		
        layout->addWidget(displayX,2,1);
                layout->addWidget(X,2,0);
        layout->addWidget(displayY,2,3);
                layout->addWidget(Y,2,2);
		
        layout->addWidget(displayJ,1,4);
                layout->addWidget(J,1,3);
		
        layout->addWidget(displayBx,3,1);
                layout->addWidget(Bx,3,0);
	layout->addWidget(displayBy,4,1);
                layout->addWidget(By,4,0);
	layout->addWidget(displayBz,5,1);
                layout->addWidget(Bz,5,0);

        this->setLayout(layout);
}

void Widget::metroButtonClicked()
{
    start_ising_window(metropolis);
}

void Widget::hbButtonClicked()
{
    start_ising_window(heatbath);
}

void Widget::wolffButtonClicked()
{
    start_ising_window(wolff);
}

void Widget::hotButtonClicked()
{
    if(hotButton->isChecked())
    {
        heat = true;
        coldButton->setChecked(false);
    }
    else
    {
        coldButton->setChecked(true);
        coldButtonClicked();
    }
}

void Widget::coldButtonClicked()
{
    if(coldButton->isChecked())
    {
        heat = false;
        hotButton->setChecked(false);
    }
    else
    {
        hotButton->setChecked(true);
        hotButtonClicked();
    }
}

void Widget::start_ising_window(SimulationType simulationType)
{
    //get values from input but catch invalid data
    float T = displayT->text().toInt(); //characters become zero
    if(isEmpty(T, *displayT, "insert finite Temp.")) return;
    
    gridSizes sizes;
    sizes.x = displayX->text().toInt(); //characters become zero
    if(isEmpty(sizes.x, *displayX, "insert grid size in x direction")) return;
    sizes.y = displayY->text().toInt(); //characters become zero
    if(isEmpty(sizes.y, *displayY, "insert grid size in y direction")) return;
    
    float J;
    coupling = displayJ->text().toFloat();
    if(isEmpty(coupling, *displayJ, "insert coupling.")) return;

    float3 B;
    B.x = displayBx->text().toFloat();
    if(isCharacter(B.x, *displayBx, "insert magn field component.")) return;
    B.y = displayBy->text().toFloat();
    if(isCharacter(B.y, *displayBy, "insert magn field component.")) return;
    B.z = displayBz->text().toFloat();
    if(isCharacter(B.z, *displayBz, "insert magn field component.")) return;

    //////////////////////////////////
    
    Window *isingWindow = new Window(sizes,
				     B, J, /*beta=*/ 1.0/(1.0*T),
			             simulationType, heat);

    isingWindow->resize(sizes.x*10, sizes.y*10);
    if(simulationType == metropolis)
      isingWindow->setWindowTitle("2D Ising Model with Metropolis algorithm");
    if(simulationType == wolff)
      isingWindow->setWindowTitle("2D Ising Model with Wolff algorithm");
    if(simulationType == heatbath)
      isingWindow->setWindowTitle("2D Ising Model with Heatbath algorithm");

    //////////////////////////////////
    isingWindow->show();

    QTimer *timer = new QTimer(this);
    connectTimerWithWindow(*timer, *isingWindow);
    timer->start(500); //update ising window every 0.5 seconds
    //////////////////////////////////
}

template<typename T>
void returnIfCharacter(T value, QLineEdit &edit, const QString errorMessage)
{
    if(!value)
    {
        if(displayB->text()=="0" || displayB->text()=="0.0")
        {
            value = 0.0f; //sets magnetic field zero if input was zero
        }
        else //error if input was character
        {
            edit.setText(errorMessage);
            return;
        }
    }
}

template<typename T>
void returnIfEmpty(T value, QLineEdit &edit, const QString errorMessage)
{
    if(!value)
    {
        displayT->setText(errorMessage);
        return;
    }
}

void connectTimerWithWindow(QTimer &timer, Window &window);
{
    connect(timer, &QTimer::timeout, window, static_cast <void (Window::*)()>(&Window::update));
    //specify correct connect function (three possibilities)
}
