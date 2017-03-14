//----------------------------------------------------------------------------------------------------------------------
/// @file main.cpp
/// @brief This is the main class that creates a GL Window
//----------------------------------------------------------------------------------------------------------------------

#include <QtGui/QApplication>
#include "MainWindow.h"

int main(int argc, char *argv[])
{
  // make an instance of the QApplication
  QApplication a(argc, argv);
  // Create a new MainWindow
  MainWindow w;
  // show it
  w.show();
  // hand control over to Qt framework
  return a.exec();
}

//----------------------------------------------------------------------------------------------------------------------
