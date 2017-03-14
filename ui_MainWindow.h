/********************************************************************************
** Form generated from reading UI file 'MainWindow.ui'
**
** Created: Mon Aug 19 18:14:41 2013
**      by: Qt User Interface Compiler version 4.8.1
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_MAINWINDOW_H
#define UI_MAINWINDOW_H

#include <QtCore/QVariant>
#include <QtGui/QAction>
#include <QtGui/QApplication>
#include <QtGui/QButtonGroup>
#include <QtGui/QCheckBox>
#include <QtGui/QComboBox>
#include <QtGui/QGridLayout>
#include <QtGui/QGroupBox>
#include <QtGui/QHeaderView>
#include <QtGui/QLabel>
#include <QtGui/QLineEdit>
#include <QtGui/QMainWindow>
#include <QtGui/QPushButton>
#include <QtGui/QSpacerItem>
#include <QtGui/QSpinBox>
#include <QtGui/QWidget>

QT_BEGIN_NAMESPACE

class Ui_MainWindow
{
public:
    QWidget *centralwidget;
    QGridLayout *s_mainGridLayout;
    QSpacerItem *horizontalSpacer;
    QGroupBox *groupBox;
    QGridLayout *gridLayout_7;
    QGroupBox *groupBox_2;
    QGridLayout *gridLayout;
    QLineEdit *m_filePath;
    QLabel *label;
    QSpacerItem *horizontalSpacer_6;
    QPushButton *m_browse;
    QGroupBox *groupBox_3;
    QGridLayout *gridLayout_3;
    QLabel *label_5;
    QComboBox *m_pointOn;
    QLabel *label_2;
    QComboBox *m_method;
    QCheckBox *m_points;
    QSpinBox *m_density;
    QLabel *label_3;
    QGroupBox *m_rayCastingGB;
    QGridLayout *gridLayout_2;
    QCheckBox *m_surfacePoints;
    QCheckBox *m_intersectionPoints;
    QGroupBox *groupBox_5;
    QGridLayout *gridLayout_5;
    QCheckBox *m_dWireframe;
    QSpacerItem *horizontalSpacer_2;
    QPushButton *m_dCompute;
    QSpacerItem *horizontalSpacer_3;
    QPushButton *m_dClear;
    QGroupBox *groupBox_6;
    QGridLayout *gridLayout_4;
    QSpacerItem *horizontalSpacer_4;
    QPushButton *m_vCompute;
    QSpacerItem *horizontalSpacer_5;
    QPushButton *m_vClear;

    void setupUi(QMainWindow *MainWindow)
    {
        if (MainWindow->objectName().isEmpty())
            MainWindow->setObjectName(QString::fromUtf8("MainWindow"));
        MainWindow->resize(1042, 723);
        centralwidget = new QWidget(MainWindow);
        centralwidget->setObjectName(QString::fromUtf8("centralwidget"));
        s_mainGridLayout = new QGridLayout(centralwidget);
        s_mainGridLayout->setObjectName(QString::fromUtf8("s_mainGridLayout"));
        horizontalSpacer = new QSpacerItem(719, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        s_mainGridLayout->addItem(horizontalSpacer, 0, 0, 1, 1);

        groupBox = new QGroupBox(centralwidget);
        groupBox->setObjectName(QString::fromUtf8("groupBox"));
        gridLayout_7 = new QGridLayout(groupBox);
        gridLayout_7->setObjectName(QString::fromUtf8("gridLayout_7"));
        groupBox_2 = new QGroupBox(groupBox);
        groupBox_2->setObjectName(QString::fromUtf8("groupBox_2"));
        gridLayout = new QGridLayout(groupBox_2);
        gridLayout->setObjectName(QString::fromUtf8("gridLayout"));
        m_filePath = new QLineEdit(groupBox_2);
        m_filePath->setObjectName(QString::fromUtf8("m_filePath"));

        gridLayout->addWidget(m_filePath, 0, 0, 1, 1);

        label = new QLabel(groupBox_2);
        label->setObjectName(QString::fromUtf8("label"));

        gridLayout->addWidget(label, 0, 1, 1, 1);

        horizontalSpacer_6 = new QSpacerItem(150, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout->addItem(horizontalSpacer_6, 1, 0, 1, 1);

        m_browse = new QPushButton(groupBox_2);
        m_browse->setObjectName(QString::fromUtf8("m_browse"));

        gridLayout->addWidget(m_browse, 1, 1, 1, 1);


        gridLayout_7->addWidget(groupBox_2, 0, 0, 1, 1);

        groupBox_3 = new QGroupBox(groupBox);
        groupBox_3->setObjectName(QString::fromUtf8("groupBox_3"));
        gridLayout_3 = new QGridLayout(groupBox_3);
        gridLayout_3->setObjectName(QString::fromUtf8("gridLayout_3"));
        label_5 = new QLabel(groupBox_3);
        label_5->setObjectName(QString::fromUtf8("label_5"));

        gridLayout_3->addWidget(label_5, 0, 0, 1, 1);

        m_pointOn = new QComboBox(groupBox_3);
        m_pointOn->setObjectName(QString::fromUtf8("m_pointOn"));

        gridLayout_3->addWidget(m_pointOn, 0, 1, 1, 2);

        label_2 = new QLabel(groupBox_3);
        label_2->setObjectName(QString::fromUtf8("label_2"));

        gridLayout_3->addWidget(label_2, 1, 0, 1, 1);

        m_method = new QComboBox(groupBox_3);
        m_method->setObjectName(QString::fromUtf8("m_method"));

        gridLayout_3->addWidget(m_method, 1, 1, 1, 2);

        m_points = new QCheckBox(groupBox_3);
        m_points->setObjectName(QString::fromUtf8("m_points"));
        m_points->setEnabled(true);
        m_points->setChecked(true);

        gridLayout_3->addWidget(m_points, 2, 0, 1, 2);

        m_density = new QSpinBox(groupBox_3);
        m_density->setObjectName(QString::fromUtf8("m_density"));
        m_density->setMaximum(1000);
        m_density->setValue(1);

        gridLayout_3->addWidget(m_density, 2, 2, 2, 1);

        label_3 = new QLabel(groupBox_3);
        label_3->setObjectName(QString::fromUtf8("label_3"));

        gridLayout_3->addWidget(label_3, 3, 0, 1, 1);

        m_rayCastingGB = new QGroupBox(groupBox_3);
        m_rayCastingGB->setObjectName(QString::fromUtf8("m_rayCastingGB"));
        m_rayCastingGB->setEnabled(false);
        gridLayout_2 = new QGridLayout(m_rayCastingGB);
        gridLayout_2->setObjectName(QString::fromUtf8("gridLayout_2"));
        m_surfacePoints = new QCheckBox(m_rayCastingGB);
        m_surfacePoints->setObjectName(QString::fromUtf8("m_surfacePoints"));
        m_surfacePoints->setEnabled(false);

        gridLayout_2->addWidget(m_surfacePoints, 0, 0, 1, 1);

        m_intersectionPoints = new QCheckBox(m_rayCastingGB);
        m_intersectionPoints->setObjectName(QString::fromUtf8("m_intersectionPoints"));
        m_intersectionPoints->setEnabled(false);

        gridLayout_2->addWidget(m_intersectionPoints, 1, 0, 1, 1);


        gridLayout_3->addWidget(m_rayCastingGB, 4, 0, 1, 3);


        gridLayout_7->addWidget(groupBox_3, 1, 0, 1, 1);

        groupBox_5 = new QGroupBox(groupBox);
        groupBox_5->setObjectName(QString::fromUtf8("groupBox_5"));
        gridLayout_5 = new QGridLayout(groupBox_5);
        gridLayout_5->setObjectName(QString::fromUtf8("gridLayout_5"));
        m_dWireframe = new QCheckBox(groupBox_5);
        m_dWireframe->setObjectName(QString::fromUtf8("m_dWireframe"));

        gridLayout_5->addWidget(m_dWireframe, 0, 0, 1, 2);

        horizontalSpacer_2 = new QSpacerItem(25, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_2, 1, 0, 1, 1);

        m_dCompute = new QPushButton(groupBox_5);
        m_dCompute->setObjectName(QString::fromUtf8("m_dCompute"));

        gridLayout_5->addWidget(m_dCompute, 1, 1, 1, 1);

        horizontalSpacer_3 = new QSpacerItem(25, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_5->addItem(horizontalSpacer_3, 1, 2, 1, 1);

        m_dClear = new QPushButton(groupBox_5);
        m_dClear->setObjectName(QString::fromUtf8("m_dClear"));

        gridLayout_5->addWidget(m_dClear, 1, 3, 1, 1);


        gridLayout_7->addWidget(groupBox_5, 2, 0, 1, 1);

        groupBox_6 = new QGroupBox(groupBox);
        groupBox_6->setObjectName(QString::fromUtf8("groupBox_6"));
        gridLayout_4 = new QGridLayout(groupBox_6);
        gridLayout_4->setObjectName(QString::fromUtf8("gridLayout_4"));
        horizontalSpacer_4 = new QSpacerItem(25, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_4, 0, 0, 1, 1);

        m_vCompute = new QPushButton(groupBox_6);
        m_vCompute->setObjectName(QString::fromUtf8("m_vCompute"));

        gridLayout_4->addWidget(m_vCompute, 0, 1, 1, 1);

        horizontalSpacer_5 = new QSpacerItem(25, 20, QSizePolicy::Expanding, QSizePolicy::Minimum);

        gridLayout_4->addItem(horizontalSpacer_5, 0, 2, 1, 1);

        m_vClear = new QPushButton(groupBox_6);
        m_vClear->setObjectName(QString::fromUtf8("m_vClear"));

        gridLayout_4->addWidget(m_vClear, 0, 3, 1, 1);


        gridLayout_7->addWidget(groupBox_6, 3, 0, 1, 1);


        s_mainGridLayout->addWidget(groupBox, 0, 1, 1, 1);

        MainWindow->setCentralWidget(centralwidget);

        retranslateUi(MainWindow);

        QMetaObject::connectSlotsByName(MainWindow);
    } // setupUi

    void retranslateUi(QMainWindow *MainWindow)
    {
        MainWindow->setWindowTitle(QApplication::translate("MainWindow", "MainWindow", 0, QApplication::UnicodeUTF8));
        groupBox->setTitle(QString());
        groupBox_2->setTitle(QApplication::translate("MainWindow", "Load Mesh :", 0, QApplication::UnicodeUTF8));
        m_filePath->setText(QApplication::translate("MainWindow", "prism.obj", 0, QApplication::UnicodeUTF8));
        label->setText(QApplication::translate("MainWindow", "(.obj)", 0, QApplication::UnicodeUTF8));
        m_browse->setText(QApplication::translate("MainWindow", "Browse", 0, QApplication::UnicodeUTF8));
        groupBox_3->setTitle(QApplication::translate("MainWindow", "Mesh Sampling", 0, QApplication::UnicodeUTF8));
        label_5->setText(QApplication::translate("MainWindow", "Point Location", 0, QApplication::UnicodeUTF8));
        m_pointOn->clear();
        m_pointOn->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "Volume", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Surface", 0, QApplication::UnicodeUTF8)
        );
        label_2->setText(QApplication::translate("MainWindow", "Method", 0, QApplication::UnicodeUTF8));
        m_method->clear();
        m_method->insertItems(0, QStringList()
         << QApplication::translate("MainWindow", "SDF", 0, QApplication::UnicodeUTF8)
         << QApplication::translate("MainWindow", "Ray Casting", 0, QApplication::UnicodeUTF8)
        );
        m_points->setText(QApplication::translate("MainWindow", "Show/Hide Points", 0, QApplication::UnicodeUTF8));
        label_3->setText(QApplication::translate("MainWindow", "Density", 0, QApplication::UnicodeUTF8));
        m_rayCastingGB->setTitle(QApplication::translate("MainWindow", "Ray Casting", 0, QApplication::UnicodeUTF8));
        m_surfacePoints->setText(QApplication::translate("MainWindow", "Surface Points", 0, QApplication::UnicodeUTF8));
        m_intersectionPoints->setText(QApplication::translate("MainWindow", "Intersection Points", 0, QApplication::UnicodeUTF8));
        groupBox_5->setTitle(QApplication::translate("MainWindow", "Delaunay", 0, QApplication::UnicodeUTF8));
        m_dWireframe->setText(QApplication::translate("MainWindow", "Wireframe", 0, QApplication::UnicodeUTF8));
        m_dCompute->setText(QApplication::translate("MainWindow", "Compute", 0, QApplication::UnicodeUTF8));
        m_dClear->setText(QApplication::translate("MainWindow", "Clear", 0, QApplication::UnicodeUTF8));
        groupBox_6->setTitle(QApplication::translate("MainWindow", "Voronoi", 0, QApplication::UnicodeUTF8));
        m_vCompute->setText(QApplication::translate("MainWindow", "Compute", 0, QApplication::UnicodeUTF8));
        m_vClear->setText(QApplication::translate("MainWindow", "Clear", 0, QApplication::UnicodeUTF8));
    } // retranslateUi

};

namespace Ui {
    class MainWindow: public Ui_MainWindow {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_MAINWINDOW_H
