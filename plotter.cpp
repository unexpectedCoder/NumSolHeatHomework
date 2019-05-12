#include "plotter.h"

#include <fstream>
#include <iostream>

#include "mainwindow.h"


using namespace QtCharts;
using namespace std;


Plotter::Plotter() : isData(false)
{
  series = new QLineSeries();
  axisX = new QValueAxis();
  axisY = new QValueAxis();
  chart = new QChart();
  chartView = new QChartView();
}


void Plotter::setData(const vector<double> &x, const vector<double> &y)
{
  if (isData)
    throw err.sendEx("data is already set");

  size_t size = x.size();
  if (size != y.size())
    throw err.sendEx("size of X != size of Y");

  for (size_t i = 0; i < size; ++i)
    series->append(x[i], y[i]);
  isData = true;
}


void Plotter::setData(const string &path, bool isHead)
{
  if (isData)
    throw err.sendEx("data is already set");

  fstream file(path.c_str(), ios_base::in);
  if (!file.is_open())
    throw err.sendEx("file is not opened");

  if (isHead)
  {
    char buf[64];
    file.getline(buf, 64, '\n');
  }

  while (true)
  {
    double x, y;
    file >> x >> y;
    if (file.eof())
      break;

    series->append(x, y);
  }
  isData = true;

  file.close();
}


void Plotter::createChart(const QString &title)
{
  if (!isData)
    throw err.sendEx("data (X and Y) is not set");

  chart->addSeries(series);
  chart->legend()->hide();
  if (!title.isEmpty())
    chart->setTitle(title);
}


void Plotter::setAxis()
{
  chart->createDefaultAxes();
}


void Plotter::setAxis(double x_min, double x_max,
                      double y_min, double y_max,
                      int x_ticks, int y_ticks)
{
  if (x_min > x_max || y_min > y_max)
    throw err.sendEx("minimum value of X or Y > maximum of it");
  if (x_ticks < 2 || y_ticks < 2)
    throw err.sendEx("ticks must be > 1");

  axisX->setMin(x_min);
  axisX->setMax(x_max);
  axisX->setTickCount(x_ticks);

  axisY->setMin(y_min);
  axisY->setMax(y_max);
  axisY->setTickCount(y_ticks);

  chart->addAxis(axisX, Qt::AlignBottom);
  chart->addAxis(axisY, Qt::AlignLeft);

  series->attachAxis(axisX);
  series->attachAxis(axisY);
}


QChartView* Plotter::getChartView()
{
  chartView->setChart(chart);
  chartView->setRenderHint(QPainter::Antialiasing);
  return chartView;
}
