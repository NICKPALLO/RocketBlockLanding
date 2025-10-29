#ifndef GRAPHIC_H
#define GRAPHIC_H

#include <QObject>

//При использовании динамической библиотеки необходимо подключить заголовочный файл,
//но не добавлять его в проект.
#include "qcustomplot.h"

inline constexpr int FIRST_GRAPH = 0;

class Graphic : public QObject
{
    Q_OBJECT   

public:
    Graphic(QCustomPlot* cPlot, uint32_t numGraps = 1);
    void AddDataToGrahp(QVector<double> x, QVector<double> y, uint32_t numGraph = 0);
    void ClearGraph(QCustomPlot* cPlot);
    void UpdateGraph(QCustomPlot* cPlot);


private:

    //Определияем указатель QCPGraph, который отвечает именно за наполнение графика.
    QList<QCPGraph*> ptrGraph;

};

#endif // GRAPHIC_H
