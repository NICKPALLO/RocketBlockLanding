#ifndef WIDGET_H
#define WIDGET_H

#include <QWidget>
#include "qcustomplot.h"
#include "landingcalc.h"
//#include <QtConcurrent>
#include "dynamiccalc.h"
#include <QFuture>
#include <QFutureWatcher>
#include "graphic.h"

QT_BEGIN_NAMESPACE
namespace Ui {
class Widget;
}
QT_END_NAMESPACE

class Widget : public QWidget
{
    Q_OBJECT

public:
    Widget(QWidget *parent = nullptr);
    ~Widget();
    LandData initializeLandDataFromUi();
    ConstData initializeDynamicDataFromUi();
public slots:
    void printLandingGraphs();
    void printDynamicGraphs();

private slots:
    void on_pushButton_clicked();

    void on_pushButton_2_clicked();

private:
    Ui::Widget *ui;
    LandingCalc landingCalc;
    DynamicCalc dynamicCalc;
    QFutureWatcher<bool> fw_land;
    QFutureWatcher<bool> fw_dynam;
    QFuture<bool> future;
    bool calculationInProgress = false;
    bool paramAvailable = false;

    void printEarth(Graphic* XY, QCustomPlot* cPlot, const QVector<double>& X, const QVector<double>& Y);
};
#endif // WIDGET_H
