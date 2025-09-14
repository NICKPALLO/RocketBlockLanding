#ifndef DYNAMICCALC_H
#define DYNAMICCALC_H

#include <QObject>
#include "dynamicRungekutt.h"
#include "aerodynamic.h"

class DynamicCalc : public QObject
{
    Q_OBJECT
public:
    explicit DynamicCalc(QObject *parent = nullptr);
    bool startCalculation(const QVector<double>& time, const QVector<double>& H_t,
                          const QVector<double>& m_t, const QVector<double>& V_t,
                          const QVector<double>& m_ok_t, const QVector<double>& m_g_t);
    void resetData();
    void setData(double h_, const ConstData& constData_);
    QVector<double>& get_y();
    QVector<double>& get_y1();
    QVector<double>& get_v();
    QVector<double>& get_v1();
    QVector<double>& get_b();
    QVector<double>& get_b1();
    QVector<double>& get_t();

signals:

private:
    ConstData constData;
    double h;
    double v;
    double v1;
    double b;
    double b1;
    double y;
    double y1;
    double t;
    QVector<double> vec_y1;
    QVector<double> vec_v1;
    QVector<double> vec_b1;
    QVector<double> vec_y;
    QVector<double> vec_v;
    QVector<double> vec_b;
    QVector<double> vec_t;
};

#endif // DYNAMICCALC_H
