#ifndef LANDINGCALC_H
#define LANDINGCALC_H
#include "matrix.h"
#include "rungekutt.h"
#include <QObject>
#include <QVector>
#include "constants.h"

#ifndef convertToMeters
#define convertToMeters *=1000
#endif

#ifndef convertToKilograms
#define convertToKilograms *=1000
#endif

#ifndef convertToNewtons
#define convertToNewtons *=-1000
#endif


class LandingCalc : public QObject
{
    Q_OBJECT
public:
    LandingCalc(QObject* parent = nullptr);
    void setParam(LandData& initialData_, float h_,double  H_atm_, double V_rq_, double delta_rq_);
    QVector<double>& getX();
    QVector<double>& getY();
    QVector<double>& getTheta();
    QVector<double>& getP_t();
    QVector<double>& getm();
    QVector<double>& getm_ok();
    QVector<double>& getm_g();
    QVector<double>& getV();
    QVector<double>& getTime();
    QVector<double>& getH();
    bool startCalculation();
    bool preCalc();
    bool braking();
    bool stopping();
    bool landing();
    void setValues();
    void modeling();
signals:
    void sendMessage(const QString& error);

private:
    LandData initialData;
    LandData currentData;
    void resetCurrentData();
    bool fuelEmpty();
    bool reached(double height);

    float h;
    double H_atm;
    double V_rq;
    double delta_rq;

    //Необходимо найти
    double V_atm;
    double t_atm;
    double theta_breaking;
    double V_stop;

    QVector<double> X_t;
    QVector<double> Y_t;
    QVector<double> theta_t;
    QVector<double> P_t;
    QVector<double> m_t;
    QVector<double> m_ok_t;
    QVector<double> m_g_t;
    QVector<double> V_t;
    QVector<double> time;
    QVector<double> H_t;

    QVector<double> K;
    QVector<double> K_t;
};

#endif // LANDINGCALC_H
