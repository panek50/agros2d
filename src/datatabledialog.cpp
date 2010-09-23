// This file is part of Agros2D.
//
// Agros2D is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 2 of the License, or
// (at your option) any later version.
//
// Agros2D is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Agros2D.  If not, see <http://www.gnu.org/licenses/>.
//
// hp-FEM group (http://hpfem.org/)
// University of Nevada, Reno (UNR) and University of West Bohemia, Pilsen
// Email: agros2d@googlegroups.com, home page: http://hpfem.org/agros2d/

#include "datatabledialog.h"

DataTableDialog::DataTableDialog(QWidget *parent) : QDialog(parent)
{
    setWindowIcon(icon(""));
    setWindowTitle(tr("Table dialog"));
    setWindowFlags(Qt::Window);

    createControls();
    highlightCurrentLineX();
    highlightCurrentLineY();

    QSettings settings;
    setMinimumSize(sizeHint());
    restoreGeometry(settings.value("DataTableDialog/Geometry", saveGeometry()).toByteArray());

    updateData();
}

DataTableDialog::~DataTableDialog()
{
    QSettings settings;
    settings.setValue("DataTableDialog/Geometry", saveGeometry());
}

void DataTableDialog::createControls()
{
    btnPlot = new QPushButton(tr("Plot data"));
    connect(btnPlot, SIGNAL(clicked()), this, SLOT(updateData()));

    txtTextX = new QTextEdit();
    connect(txtTextX, SIGNAL(cursorPositionChanged()), this, SLOT(highlightCurrentLineX()));
    txtTextY = new QTextEdit();
    connect(txtTextY, SIGNAL(cursorPositionChanged()), this, SLOT(highlightCurrentLineY()));

    QGridLayout *textLayout = new QGridLayout();
    textLayout->addWidget(new QLabel("X"), 0, 0);
    textLayout->addWidget(new QLabel("Y"), 0, 1);
    textLayout->addWidget(txtTextX, 1, 0);
    textLayout->addWidget(txtTextY, 1, 1);
    textLayout->addWidget(btnPlot, 2, 0, 1, 2);

    QGroupBox *tableGroup = new QGroupBox(tr("Table values"));
    // tableWidget->setLayout(tableLayout);
    tableGroup->setLayout(textLayout);
    tableGroup->setMinimumWidth(250);
    tableGroup->setMaximumWidth(250);

    // chart
    chartValue = new Chart(this);

    // curve
    curveValue = new QwtPlotCurve();
    curveValue->setRenderHint(QwtPlotItem::RenderAntialiased);
    curveValue->setPen(QPen(Qt::blue));
    curveValue->setCurveAttribute(QwtPlotCurve::Inverted);
    curveValue->setYAxis(QwtPlot::yLeft);
    curveValue->setTitle(tr("value"));
    curveValue->attach(chartValue);

    // labels
    QwtText textErrorBottom(tr("x"));
    textErrorBottom.setFont(QFont("Helvetica", 10, QFont::Normal));
    chartValue->setAxisTitle(QwtPlot::xBottom, textErrorBottom);

    QwtText textErrorLeft(tr("y"));
    textErrorLeft.setFont(QFont("Helvetica", 10, QFont::Normal));
    chartValue->setAxisTitle(QwtPlot::yLeft, textErrorLeft);

    chartDerivative = new Chart(this);

    // curve
    curveDerivative = new QwtPlotCurve();
    curveDerivative->setRenderHint(QwtPlotItem::RenderAntialiased);
    curveDerivative->setPen(QPen(Qt::blue));
    curveDerivative->setCurveAttribute(QwtPlotCurve::Inverted);
    curveDerivative->setYAxis(QwtPlot::yLeft);
    curveDerivative->setTitle(tr("derivative"));
    curveDerivative->attach(chartDerivative);

    // labels
    QwtText textDerivativeBottom(tr("x"));
    textDerivativeBottom.setFont(QFont("Helvetica", 10, QFont::Normal));
    chartDerivative->setAxisTitle(QwtPlot::xBottom, textDerivativeBottom);

    QwtText textDerivativeLeft(tr("dy/dx"));
    textDerivativeLeft.setFont(QFont("Helvetica", 10, QFont::Normal));
    chartDerivative->setAxisTitle(QwtPlot::yLeft, textDerivativeLeft);

    QVBoxLayout *chartLayout = new QVBoxLayout();
    chartLayout->addWidget(chartValue);
    chartLayout->addWidget(chartDerivative);

    // main layout
    QHBoxLayout *mainLayout = new QHBoxLayout();
    mainLayout->addWidget(tableGroup);
    mainLayout->addLayout(chartLayout);

    // dialog buttons
    QDialogButtonBox *buttonBox = new QDialogButtonBox(QDialogButtonBox::Ok | QDialogButtonBox::Cancel);
    connect(buttonBox, SIGNAL(accepted()), this, SLOT(doAccept()));
    connect(buttonBox, SIGNAL(rejected()), this, SLOT(doReject()));

    QVBoxLayout *layout = new QVBoxLayout();
    layout->addLayout(mainLayout);
    layout->addWidget(buttonBox);

    setLayout(layout);
}

void DataTableDialog::updateData()
{
    QStringList listX = txtTextX->toPlainText().trimmed().
                        replace(" ", "").
                        replace("\n", ";").
                        split(";", QString::KeepEmptyParts);
    QStringList listY = txtTextY->toPlainText().trimmed().
                        replace(" ", "").
                        replace("\n", ";").
                        split(";", QString::KeepEmptyParts);

    if (listX.count() != listY.count())
    {
        QMessageBox::critical(QApplication::activeWindow(), tr("Table"), tr("Length of both datasets must be same."));
        return;
    }

    dataTable.clear();
    for (int i = 0; i < listX.count(); i++)
    {
        bool okX;
        bool okY;

        double x = QString(listX[i]).toDouble(&okX);
        double y = QString(listY[i]).toDouble(&okY);

        if (okX && okY)
            dataTable.add(x, y);
    }

    // reload and sort data
    txtTextX->clear();
    txtTextY->clear();

    /*
    DataTableRow *data = dataTable.data();
    while (data)
    {
        txtTextX->append(QString("%1").arg(data->key));
        txtTextY->append(QString("%1").arg(data->value));

        // next row
        data = data->next;
    }
    */
    gotoLineX(0);
    gotoLineY(0);

    // plot chart
    int count = 100;
    double min = dataTable.min_key();
    double max = dataTable.max_key();

    double *xval = new double[count];
    double *yvalValue = new double[count];
    double *yvalDerivative = new double[count];

    for (int i = 0; i < count; i++)
    {
        double key = min + i * (max - min) / (count+1);

        xval[i] = key;
        yvalValue[i] = dataTable.value(key);
        yvalDerivative[i] = dataTable.derivative(key);
    }

    // plot value
    bool doReplotValue = chartValue->autoReplot();
    chartValue->setAutoReplot(false);

    curveValue->setData(xval, yvalValue, count);

    chartValue->setAutoReplot(doReplotValue);
    chartValue->replot();

    // plot derivative
    bool doReplotDerivative = chartValue->autoReplot();
    chartDerivative->setAutoReplot(false);

    curveDerivative->setData(xval, yvalDerivative, count);

    chartDerivative->setAutoReplot(doReplotDerivative);
    chartDerivative->replot();

    delete[] xval;
    delete[] yvalValue;
    delete[] yvalDerivative;
}

void DataTableDialog::doAccept()
{
    // if (save()) accept();
}

void DataTableDialog::doReject()
{
    reject();
}

void DataTableDialog::highlightCurrentLineX()
{
    QList<QTextEdit::ExtraSelection> selections;

    QTextEdit::ExtraSelection selection;
    QColor lineColor = QColor(Qt::yellow).lighter(180);

    selection.format.setBackground(lineColor);
    selection.format.setProperty(QTextFormat::FullWidthSelection, true);
    selection.cursor = txtTextX->textCursor();
    selection.cursor.clearSelection();
    selections.append(selection);

    txtTextX->setExtraSelections(selections);
    gotoLineY(txtTextX->textCursor().blockNumber() + 1);
}

void DataTableDialog::highlightCurrentLineY()
{
    QList<QTextEdit::ExtraSelection> selections;

    QTextEdit::ExtraSelection selection;
    QColor lineColor = QColor(Qt::yellow).lighter(180);

    selection.format.setBackground(lineColor);
    selection.format.setProperty(QTextFormat::FullWidthSelection, true);
    selection.cursor = txtTextY->textCursor();
    selection.cursor.clearSelection();
    selections.append(selection);

    txtTextY->setExtraSelections(selections);
    // gotoLineX(txtTextY->textCursor().blockNumber());
}

void DataTableDialog::gotoLineX(int line)
{
    if (line >= 0 && line <= txtTextX->document()->blockCount())
    {
        int pos = txtTextX->document()->findBlockByNumber(line - 1).position();
        QTextCursor cur = txtTextX->textCursor();
        cur.setPosition(pos, QTextCursor::MoveAnchor);
        txtTextX->setTextCursor(cur);
        txtTextX->ensureCursorVisible();
        highlightCurrentLineX();
    }
}

void DataTableDialog::gotoLineY(int line)
{
    if (line >= 0 && line <= txtTextY->document()->blockCount())
    {
        int pos = txtTextY->document()->findBlockByNumber(line - 1).position();
        QTextCursor cur = txtTextY->textCursor();
        cur.setPosition(pos, QTextCursor::MoveAnchor);
        txtTextY->setTextCursor(cur);
        txtTextY->ensureCursorVisible();
        highlightCurrentLineY();
    }
}
