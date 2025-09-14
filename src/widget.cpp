#include "widget.h"
#include <QtConcurrent>
#include "qcustomplot.h"
#include "../form/ui_widget.h"
#include "rungekutt.h"
#include"dynamicRungekutt.h"
#include <algorithm>

Widget::Widget(QWidget *parent)
    : QWidget(parent)
    , ui(new Ui::Widget)
{
    ui->setupUi(this);
    QObject::connect(&landingCalc, &LandingCalc::sendMessage, this, [&](QString message){
        ui->lab_inf->setText(message);
    });
    connect(&fw_land,&QFutureWatcher<bool>::finished,this,&Widget::printLandingGraphs);
    connect(&fw_dynam,&QFutureWatcher<bool>::finished,this,&Widget::printDynamicGraphs);
}

Widget::~Widget()
{
    delete ui;
}

LandData Widget::initializeLandDataFromUi()
{
    double Xc = ui->le_Xc->text().toDouble();
    double Yc = ui->le_Yc->text().toDouble();
    double Zc = ui->le_Zc->text().toDouble();
    double Vx_c = ui->le_Vx_c->text().toDouble();
    double Vy_c = ui->le_Vy_c->text().toDouble();
    double Vz_c = ui->le_Vz_c->text().toDouble();
    double m = ui->le_m->text().toDouble();
    double theta = ui->le_theta->text().toDouble();
    double psi = ui->le_psi->text().toDouble();
    double gamma = ui->le_gamma->text().toDouble();
    double Az = ui->le_Az->text().toDouble();
    double ga = ui->le_ga->text().toDouble();
    double la = ui->le_la->text().toDouble();
    double a = ui->le_a->text().toDouble();
    double B = ui->le_B->text().toDouble();
    double Sm = ui->le_Sm->text().toDouble();
    double mt_ok = ui->le_mt_ok->text().toDouble();
    double dm_ok = ui->le_dm_ok->text().toDouble();
    double mt_g = ui->le_mt_g->text().toDouble();
    double dm_g = ui->le_dm_g->text().toDouble();
    double P = ui->le_P->text().toDouble();
    Xc convertToMeters;
    Yc convertToMeters;
    Zc convertToMeters;
    m convertToKilograms;
    mt_ok convertToKilograms;
    dm_ok convertToKilograms;
    mt_g convertToKilograms;
    dm_g convertToKilograms;
    P convertToNewtons;
    double nn = 90 - Az;
    vec V_c(Vx_c, Vy_c, Vz_c);
    vec r_c(Xc, Yc, Zc);
    vec r_m = get_rm_from_rc(r_c, nn, la);

    LandData initialData;
    initialData.t = 0;
    initialData.m = m;
    initialData.gamma = gamma;
    initialData.psi = psi;
    initialData.theta = theta;
    initialData.a = a;
    initialData.B = B;
    initialData.r_c = r_c;
    initialData.r_m = r_m;
    initialData.V_c = V_c;
    initialData.mt_ok = mt_ok;
    initialData.mt_g = mt_g;
    initialData.kp = 1;
    initialData.P = P;
    initialData.Sm = Sm;
    initialData.nn = nn;
    initialData.la = la;
    initialData.mt_ok_n = mt_ok;
    initialData.dm_ok = dm_ok;
    initialData.mt_g_n = mt_g;
    initialData.dm_g = dm_g;

    initialData.tpn[0] = 1;
    initialData.tpn[1] = 1;
    initialData.tpn[2] = 1;
    initialData.tpend[0] = -1;
    initialData.tpend[1] = -1;
    initialData.tpend[2] = -1;

    return initialData;
}

ConstData Widget::initializeDynamicDataFromUi()
{
    ConstData initialData;
    initialData.a0 = ui->le_a0->text().toDouble();
    initialData.a1 = ui->le_a1->text().toDouble();
    initialData.a2 = ui->le_a2->text().toDouble();
    initialData.a3 = ui->le_a3->text().toDouble();
    initialData.t1 = ui->le_t1->text().toDouble();
    initialData.t2 = ui->le_t2->text().toDouble();
    //initialData.P = ui->le_P->text().toDouble();
    initialData.P = 98060;
    initialData.Sm = ui->le_Sm->text().toDouble();
    initialData.D = ui->le_Drb->text().toDouble();
    initialData.L = ui->le_Lrb->text().toDouble();
    initialData.mk = ui->le_mk->text().toDouble();
    initialData.p_ok = ui->le_pok->text().toDouble();
    initialData.p_g = ui->le_pg->text().toDouble();
    initialData.Xk = ui->le_Xk->text().toDouble();
    initialData.Xdn_ok = ui->le_Xdok->text().toDouble();
    initialData.Xdn_g = ui->le_Xdg->text().toDouble();
    return initialData;
}

void Widget::printLandingGraphs()
{
    bool res = future.result();
    if(res)
    {
        calculationInProgress = false;
        Graphic grafic_XY(ui->XY_graph,2);
        Graphic grafic_Theta(ui->theta_graph);
        Graphic grafic_P(ui->P_graph);
        Graphic grafic_m(ui->m_graph);
        Graphic grafic_m_ok(ui->m_ok_graph);
        Graphic grafic_m_g(ui->m_g_graph);
        Graphic grafic_V(ui->V_graph);

        grafic_XY.ClearGraph(ui->XY_graph);
        grafic_Theta.ClearGraph(ui->theta_graph);
        grafic_P.ClearGraph(ui->P_graph);
        grafic_m.ClearGraph(ui->m_graph);
        grafic_m_ok.ClearGraph(ui->m_ok_graph);
        grafic_m_g.ClearGraph(ui->m_g_graph);
        grafic_V.ClearGraph(ui->V_graph);


        // int maxX = *(std::max_element(landingCalc.getX().begin(),landingCalc.getX().end()));
        // int minX = *(std::min_element(landingCalc.getX().begin(),landingCalc.getX().end()));
        // int maxY = *(std::max_element(landingCalc.getY().begin(),landingCalc.getY().end()));
        // int minY = *(std::min_element(landingCalc.getY().begin(),landingCalc.getY().end()));

        // int max = maxX>maxY ? maxX : maxY;
        // int min = minY<minX ? minY : minX;
        // max*=1.3;
        // min = min>0 ? -10 : min*1.3;

        // QVector<double> XR;
        // QVector<double> YR;

        // double Xr = min;
        // double Yr;
        // while(Xr<max)
        // {
        //     Yr = sqrt((6371.1*6371.1)-Xr*Xr)-6371.1;
        //     XR.push_back(Xr);
        //     YR.push_back(Yr);
        //     Xr+=10;
        // }


        // ui->XY_graph->xAxis->setRange(min,max);
        // ui->XY_graph->yAxis->setRange(min,max);
        grafic_XY.AddDataToGrahp(landingCalc.getX(),landingCalc.getY(),0);
        // grafic_XY.AddDataToGrahp(XR,YR,1);


        grafic_Theta.AddDataToGrahp(landingCalc.getTime(),landingCalc.getTheta());
        grafic_P.AddDataToGrahp(landingCalc.getTime(),landingCalc.getP_t());
        grafic_m.AddDataToGrahp(landingCalc.getTime(),landingCalc.getm());
        grafic_m_ok.AddDataToGrahp(landingCalc.getTime(),landingCalc.getm_ok());
        grafic_m_g.AddDataToGrahp(landingCalc.getTime(),landingCalc.getm_g());
        grafic_V.AddDataToGrahp(landingCalc.getTime(),landingCalc.getV());

        //grafic_XY.UpdateGraph(ui->XY_graph);
        printEarth(&(grafic_XY),ui->XY_graph,landingCalc.getX(), landingCalc.getY());
        grafic_Theta.UpdateGraph(ui->theta_graph);
        grafic_P.UpdateGraph(ui->P_graph);
        grafic_m.UpdateGraph(ui->m_graph);
        grafic_m_ok.UpdateGraph(ui->m_ok_graph);
        grafic_m_g.UpdateGraph(ui->m_g_graph);
        grafic_V.UpdateGraph(ui->V_graph);

        calculationInProgress = false;
        paramAvailable = true;
    }
    else
    {
        calculationInProgress = false;
        paramAvailable = false;
    }
}

void Widget::printDynamicGraphs()
{
    bool res = future.result();
    if(res)
    {
        Graphic grafic_y(ui->yt_graph);
        Graphic grafic_y1(ui->y1t_graph);
        Graphic grafic_v(ui->vt_graph);
        Graphic grafic_v1(ui->v1t_graph);
        Graphic grafic_b(ui->bt_graph);
        Graphic grafic_b1(ui->b1t_graph);

        grafic_y.ClearGraph(ui->yt_graph);
        grafic_y1.ClearGraph(ui->y1t_graph);
        grafic_v.ClearGraph(ui->vt_graph);
        grafic_v1.ClearGraph(ui->v1t_graph);
        grafic_b.ClearGraph(ui->bt_graph);
        grafic_b1.ClearGraph(ui->b1t_graph);

        grafic_y.AddDataToGrahp(dynamicCalc.get_t(),dynamicCalc.get_y());
        grafic_y1.AddDataToGrahp(dynamicCalc.get_t(),dynamicCalc.get_y1());
        grafic_v.AddDataToGrahp(dynamicCalc.get_t(),dynamicCalc.get_v());
        grafic_v1.AddDataToGrahp(dynamicCalc.get_t(),dynamicCalc.get_v1());
        grafic_b.AddDataToGrahp(dynamicCalc.get_t(),dynamicCalc.get_b());
        grafic_b1.AddDataToGrahp(dynamicCalc.get_t(),dynamicCalc.get_b1());

        grafic_y.UpdateGraph(ui->yt_graph);
        grafic_y1.UpdateGraph(ui->y1t_graph);
        grafic_v.UpdateGraph(ui->vt_graph);
        grafic_v1.UpdateGraph(ui->v1t_graph);
        grafic_b.UpdateGraph(ui->bt_graph);
        grafic_b1.UpdateGraph(ui->b1t_graph);
    }
    calculationInProgress = false;
}

void Widget::on_pushButton_clicked()
{
    if(!calculationInProgress)
    {
        calculationInProgress = true;
        LandData initialData = initializeLandDataFromUi();
        float h = ui->le_h->text().toDouble();
        double H_atm = ui->le_H_atm->text().toDouble();
        H_atm *=1000;
        double V_rq = ui->le_V_rq->text().toDouble();
        double delta_rq = ui->le_delta_rq->text().toDouble();
        landingCalc.setParam(initialData, h, H_atm, V_rq, delta_rq);
        future = QtConcurrent::run(&LandingCalc::startCalculation,&landingCalc);
        fw_land.setFuture(future);
    }
}


void Widget::on_pushButton_2_clicked()
{
    if(!calculationInProgress && paramAvailable)
    {
        calculationInProgress = true;
        ConstData constData = initializeDynamicDataFromUi();
        double h =ui->le_h_dynam->text().toDouble();
        dynamicCalc.setData(h,constData);
        future = QtConcurrent::run(&DynamicCalc::startCalculation,&dynamicCalc,std::cref(landingCalc.getTime()),
                                   std::cref(landingCalc.getH()),std::cref(landingCalc.getm()),std::cref(landingCalc.getV()),
                                   std::cref(landingCalc.getm_ok()),std::cref(landingCalc.getm_g()));
        fw_dynam.setFuture(future);
    }
}

void Widget::printEarth(Graphic* XY, QCustomPlot* cPlot,const QVector<double>& X, const QVector<double>& Y)
{
    int maxX = *(std::max_element(X.begin(),X.end()));
    int minX = *(std::min_element(X.begin(),X.end()));
    int maxY = *(std::max_element(Y.begin(),Y.end()));
    int minY = *(std::min_element(Y.begin(),Y.end()));

    int max = maxX>maxY ? maxX : maxY;
    int min = minY<minX ? minY : minX;
    max*=1.3;
    min = min>0 ? -10 : min*1.3;

    QVector<double> XR;
    QVector<double> YR;

    double Xr = min;
    double Yr;
    while(Xr<max)
    {
        Yr = sqrt((6371.1*6371.1)-Xr*Xr)-6371.1;
        XR.push_back(Xr);
        YR.push_back(Yr);
        Xr+=10;
    }
    XY->AddDataToGrahp(XR,YR,1);
    cPlot->xAxis->setRange(min,max);
    cPlot->yAxis->setRange(min,max);
}

