#ifndef MAINWINDOW_H
#define MAINWINDOW_H

//----------------------------------------------------------------------------------------------------------------------
/// @file MainWindow.h
/// @brief This is the MainWindow Class which is generated by the Ui file, if we wish to add anything to the main
/// Ui we add it here
/// @author Maria Vineeta Bagya Seelan
/// @version 1.0
/// @date 16/08/13
/// Revision History :
/// Initial Version 16/08/13
/// @class MainWindow
/// @brief the main re-sizable window which contains a GLWindow widget used to hold our
/// basic gl applications
//----------------------------------------------------------------------------------------------------------------------

#include "Renderer.h"
#include <QtGui/QMainWindow>

/// @namespace Ui our Ui namespace created from the MainWindow classMain
namespace Ui {
              class MainWindow;
             }

class MainWindow : public QMainWindow
{
    Q_OBJECT

protected :
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief override the keyPressEvent inherited from QObject so we can handle key presses.
  /// @param [in] _event the event to process
  //----------------------------------------------------------------------------------------------------------------------
  void keyPressEvent(
                     QKeyEvent *_event
                    );

public:
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief constructor
  /// @param _parent the parent window the for this window
  //----------------------------------------------------------------------------------------------------------------------
  MainWindow(
              QWidget *_parent = 0
            );
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief  dtor free up the GLWindow and all resources
  //----------------------------------------------------------------------------------------------------------------------
  ~MainWindow();

private slots :
  /// @brief A slot to browse for the Obj
  /// Button
  void loadObjFile();
  /// @brief A slot to update the method for the obj
  /// Combo Box
  void methodChanged( int _method );
  /// @brief A slot to update the point location Surface/Volume
  /// Combo Box
  void setPointType( int _type );

private:
  //----------------------------------------------------------------------------------------------------------------------
  ///  @brief our gl window created in GLWindow.h
  //----------------------------------------------------------------------------------------------------------------------
  Renderer *m_gl;
  //----------------------------------------------------------------------------------------------------------------------
  /// @brief Our MainWindow UI
  //----------------------------------------------------------------------------------------------------------------------
  Ui::MainWindow *m_ui;

};

#endif // MAINWINDOW_H
//------------------------------------------------------------------------------------------------------------------------
