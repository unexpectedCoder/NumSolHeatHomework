#ifndef PLOTTER_H
#define PLOTTER_H

#include <QWidget>

#include <QtCharts/QChart>
#include <QtCharts/QChartView>
#include <QtCharts/QLineSeries>
#include <QtCharts/QValueAxis>

#include <QPainter>
#include <QPen>
#include <QFont>
#include <QString>

#include <vector>

#include "err.h"


class Plotter
{
private:
  Error err;

  QtCharts::QLineSeries *series;
  QtCharts::QValueAxis *axisX, *axisY;
  QtCharts::QChart *chart;
  QtCharts::QChartView *chartView;

  bool isData;

public:
  Plotter();

  void setData(const std::vector<double> &x, const std::vector<double> &y);
  void setData(const std::string &path, bool isHead);
  void createChart(const QString &title = "");
  void setAxis();
  void setAxis(double x_min, double x_max,
               double y_min, double y_max,
               int x_ticks, int y_ticks);
  QtCharts::QChartView* getChartView();
};


#endif // PLOTTER_H
