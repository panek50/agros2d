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

#ifndef ACOUSTIC_H
#define ACOUSTIC_H

#include "util.h"
#include "module.h"

struct ModuleAcoustic : public Hermes::Module::ModuleAgros
{
    Q_OBJECT
public:
    ModuleAcoustic(ProblemType problemType, AnalysisType analysisType) : Hermes::Module::ModuleAgros(problemType, analysisType) {}

    inline int number_of_solution() const { return 2; }
    bool has_nonlinearity() const { return false; }

    LocalPointValue *local_point_value(const Point &point);
    SurfaceIntegralValue *surface_integral_value();
    VolumeIntegralValue *volume_integral_value();

    ViewScalarFilter *view_scalar_filter(Hermes::Module::LocalVariable *physicFieldVariable,
                                         PhysicFieldVariableComp physicFieldVariableComp);

    Hermes::vector<SolutionArray *> solve(ProgressItemSolve *progressItemSolve);
    void update_time_functions(double time);

    // rewrite
    void readBoundaryFromDomElement(QDomElement *element);
    void writeBoundaryToDomElement(QDomElement *element, SceneBoundary *marker);
    void readMaterialFromDomElement(QDomElement *element);
    void writeMaterialToDomElement(QDomElement *element, SceneMaterial *marker);

    SceneBoundary *newBoundary();
    SceneBoundary *newBoundary(PyObject *self, PyObject *args);
    SceneBoundary *modifyBoundary(PyObject *self, PyObject *args);
    SceneMaterial *newMaterial();
    SceneMaterial *newMaterial(PyObject *self, PyObject *args);
    SceneMaterial *modifyMaterial(PyObject *self, PyObject *args);
};

// *******************************************************************************************

class ParserAcoustic : public Parser
{
public:
    double prho;
    double pv;

    void setParserVariables(SceneMaterial *material);
};

class LocalPointValueAcoustic : public LocalPointValue
{
public:
    LocalPointValueAcoustic(const Point &point);
};

class SurfaceIntegralValueAcoustic : public SurfaceIntegralValue
{
public:
    SurfaceIntegralValueAcoustic();
};

class VolumeIntegralValueAcoustic : public VolumeIntegralValue
{
public:
    VolumeIntegralValueAcoustic();
};

class ViewScalarFilterAcoustic : public ViewScalarFilter
{
public:
    ViewScalarFilterAcoustic(Hermes::vector<MeshFunction *> sln,
                                   std::string expression);
};

// *******************************************************************************************

class SceneBoundaryAcoustic : public SceneBoundary
{
public:
    Value value_real;

    SceneBoundaryAcoustic(const QString &name, PhysicFieldBC type, Value value_real);
    SceneBoundaryAcoustic(const QString &name, PhysicFieldBC type);

    QString script();
    QMap<QString, QString> data();
    int showDialog(QWidget *parent);
};

class SceneMaterialAcoustic : public SceneMaterial
{
public:
    Value density;
    Value speed;

    SceneMaterialAcoustic(const QString &name, Value density, Value speed);

    QString script();
    QMap<QString, QString> data();
    int showDialog(QWidget *parent);
};

class SceneBoundaryAcousticDialog : public SceneBoundaryDialog
{
    Q_OBJECT

public:
    SceneBoundaryAcousticDialog(SceneBoundaryAcoustic *boundary, QWidget *parent);

protected:
    void createContent();

    void load();
    bool save();

private:
    QLabel *lblValueUnit;
    QComboBox *cmbType;
    ValueLineEdit *txtValue;

private slots:
    void doTypeChanged(int index);
};

class SceneMaterialAcousticDialog : public SceneMaterialDialog
{
    Q_OBJECT

public:
    SceneMaterialAcousticDialog(SceneMaterialAcoustic *material, QWidget *parent);

protected:
    void createContent();

    void load();
    bool save();

private:
    ValueLineEdit *txtDensity;
    ValueLineEdit *txtSpeed;
};

#endif // ACOUSTIC_H
