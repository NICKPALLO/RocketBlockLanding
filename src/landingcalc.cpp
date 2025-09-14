#include "landingcalc.h"

LandingCalc::LandingCalc(QObject* parent) : QObject(parent) {}

void LandingCalc::setParam(LandData &initialData_, float h_, double H_atm_, double V_rq_, double delta_rq_)
{
    initialData = initialData_;
    currentData = initialData_;
    h = h_;
    H_atm = H_atm_;
    V_rq = V_rq_;
    delta_rq = delta_rq_;
    X_t.clear();
    Y_t.clear();
    theta_t.clear();
    P_t.clear();
    m_t.clear();
    m_ok_t.clear();
    m_g_t.clear();
    V_t.clear();
    time.clear();
    H_t.clear();
    K.clear();
    K_t.clear();
}

QVector<double>& LandingCalc::getX()
{
    return X_t;
}

QVector<double>& LandingCalc::getY()
{
    return Y_t;
}

QVector<double> &LandingCalc::getTheta()
{
    return theta_t;
}

QVector<double> &LandingCalc::getP_t()
{
    return P_t;
}

QVector<double> &LandingCalc::getm()
{
    return m_t;
}

QVector<double> &LandingCalc::getm_ok()
{
    return m_ok_t;
}

QVector<double> &LandingCalc::getm_g()
{
    return m_g_t;
}

QVector<double> &LandingCalc::getV()
{
    return V_t;
}

QVector<double> &LandingCalc::getTime()
{
    return time;
}

QVector<double> &LandingCalc::getH()
{
    return H_t;
}

void LandingCalc::resetCurrentData()
{
    std::array<double,AmountPulses> tpnTEMP = std::move(currentData.tpn);
    std::array<double,AmountPulses> tpendTEMP = std::move(currentData.tpend);
    currentData = initialData;
    currentData.tpn = std::move(tpnTEMP);
    currentData.tpend = std::move(tpendTEMP);
}

bool LandingCalc::fuelEmpty()
{
    if(currentData.mt_g>=0 && currentData.mt_ok>=0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

bool LandingCalc::reached(double height)
{
    if(currentData.r_m.size() < Rz+height)
    {
        return true;
    }
    else
    {
        return false;
    }
}

bool LandingCalc::preCalc()
{
    int max = 0;
    int min = 0;
    bool first = true;
    do
    {
        rungeKutt(currentData, h);
        //для рассчета тормозного импульса (Breaking)
        if(reached(H_atm) && currentData.V_c.y<=0 && first)
        {
            V_atm = currentData.V_c.size();
            t_atm = currentData.t;
            first = false;
        }


    }while(!reached(0));


    return true;
}

bool LandingCalc::braking()
{
    resetCurrentData();
    currentData.tpn[0] = t_atm;
    currentData.tpend[0] = 100000; //начнется, но никогда не закончится
    bool braking = false;
    bool first = true;
    theta_breaking = 0;
    if(V_atm>V_rq);
    {
        while(!braking && currentData.tpn[0]>=0 && !fuelEmpty())
        {
            resetCurrentData();

            currentData.tpn[0] -= h;
            first = true;
            while((!reached(H_atm) || currentData.V_c.y>0) && !fuelEmpty())
            {
                rungeKutt(currentData, h);
                if (currentData.t>currentData.tpn[0] && first)
                {
                    theta_breaking = getThetafromVc(currentData.V_c);
                    theta_breaking = 180-abs(theta_breaking);
                    first = false;
                }
            }

            if(currentData.V_c.size()<V_rq)
            {
                braking = true;
                currentData.tpend[0] = currentData.t;
                while(!reached(10))
                {
                    rungeKutt(currentData, h);
                }
                currentData.tpn[1] = currentData.t;
                currentData.tpend[1] = 100000;
                V_stop = currentData.V_c.size();
            }
        }
    }
    return braking;
}

bool LandingCalc::stopping()
{
    resetCurrentData();

    double VertV = 0;
    bool stop = false;
    bool flyUp = false;
    if(!reached(10) && V_stop>0)
    {
        while(!stop && currentData.tpn[1]>=currentData.tpend[0] && !fuelEmpty())
        {
            resetCurrentData();
            currentData.tpn[1] -= h;
            flyUp = false;
            while(!reached(5) && !fuelEmpty())
            {
                rungeKutt(currentData, h);
                VertV=getVertV(currentData);
                if (VertV>0 && currentData.t>currentData.tpn[1])
                {
                    flyUp = true;
                    break;
                }

            }
            if((abs(VertV)<5 && reached(15) && !reached(0)) || flyUp)
            {
                stop = true;
                currentData.tpend[1] = currentData.t;
            }
        }
    }
    return stop;
}

bool LandingCalc::landing()
{
    bool landing = false;

    resetCurrentData();

    while (!reached(0))
    {
        rungeKutt(currentData, h);
        if(currentData.t >= currentData.tpend[1] && currentData.V_c.y<0)
        {
            break;
        }
    }
    currentData.tpn[2] = currentData.t;
    currentData.tpend[2]=currentData.t+1000;
    LandData previousData = currentData;
    bool skipStep = false;
    bool returnToPrevStep = false;
    double P_upr=-currentData.m * 9.81;
    double VertV = 0;
    while(!reached(0) && !fuelEmpty())
    {
        if(abs(P_upr)<abs(currentData.P))
        {
            currentData.kp=P_upr/currentData.P;
        }
        else
        {
            currentData.kp = 1;
        }
        rungeKutt(currentData, h);
        VertV = getVertV(currentData);
        if((VertV>=-2 && VertV<-1) || skipStep)
        {
            previousData = currentData;
            K.push_back(currentData.kp);
            K_t.push_back(currentData.t-h);
        }
        else
        {
            returnToPrevStep=true;
        }
        if(!(VertV>=-2 && VertV<-1))
        {
            if(VertV>=-1)
            {
                if (P_upr<currentData.P * 0.001)
                {
                    P_upr -= currentData.P * 0.001;
                    skipStep = false;
                }
                else
                {
                    P_upr=0;
                    skipStep = true;
                }
            }
            else
            {
                if (P_upr>currentData.P*0.999)
                {
                    P_upr += currentData.P * 0.001;
                    skipStep = false;
                }
                else
                {
                    P_upr=currentData.P;
                    skipStep = true;
                }
            }
        }
        if(returnToPrevStep)
        {
            currentData = previousData;
            returnToPrevStep=false;
        }
    }
    if(!fuelEmpty())
    {
        landing = true;
        currentData.tpend[2] = currentData.t;
    }

    return landing;
}

void LandingCalc::setValues()
{
    X_t.push_back(currentData.r_c.x/1000);
    Y_t.push_back(currentData.r_c.y/1000);
    P_t.push_back(-(getP(currentData.t,currentData).x));
    time.push_back(currentData.t);
    m_t.push_back(currentData.m);
    m_ok_t.push_back(currentData.mt_ok);
    m_g_t.push_back(currentData.mt_g);
    V_t.push_back(sqrt(currentData.V_c.x*currentData.V_c.x+currentData.V_c.y*currentData.V_c.y));
    H_t.push_back((currentData.r_m.size()-Rz)/1000);
    theta_t.push_back(getRealTheta(theta_breaking,initialData.theta,currentData));
}

void LandingCalc::modeling()
{
    resetCurrentData();
    do
    {
        setValues();
        for(int i = 0; i <K_t.size();++i)
        {
            if(currentData.t == K_t[i])
            {
                currentData.kp = K[i];
            }
        }
        rungeKutt(currentData, h);
    }while(!reached(0));
}


bool LandingCalc::startCalculation()
{
    bool land = false;
    bool stop = false;
    emit sendMessage("Идет предварительный рассчет...");
    preCalc();
    emit sendMessage("Идет рассчет торможения в плотных слоях атмосферы...");
    bool brake = braking();
    if(brake)
    {
        emit sendMessage("Идет рассчет процесса остановки...");
        stop = stopping();
        if(stop)
        {
            emit sendMessage("Идет рассчет процесса вертикальной посадки...");
            land = landing();
            if(land)
            {
                modeling();
            }
        }
    }
    if(!brake)
    {
        if(fuelEmpty())
        {
            emit sendMessage("не хватило топлива на этапе торможения в плотных слоях атмосферы");
        }
        else
        {
            emit sendMessage("не хватило времени для выдачи импульса на этапе торможения в плотных слоях атмосферы");
        }
    }
    else if(!stop)
    {
        if(fuelEmpty())
        {
            emit sendMessage("не хватило топлива на этапе остановки ракетного блока");
        }
        else
        {
            emit sendMessage("не хватило времени для выдачи импульса на этапе остановки ракетного блока");
        }
    }
    else if(!land)
    {
        if(fuelEmpty())
        {
            emit sendMessage("не хватило топлива на этапе вертикального спуска");
        }
        else
        {
            emit sendMessage("не хватило времени для выдачи импульса на этапе вертикального спуска");
        }
    }
    else
    {
        emit sendMessage("Рассчет успешно завершен");
    }
    return land;
}

