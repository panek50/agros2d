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

#include "block.h"
#include "util/global.h"

#include "field.h"
#include "problem.h"
#include "coupling.h"
#include "scene.h"
#include "logview.h"
#include "solver.h"
#include "module.h"
#include "module_agros.h"

Block::Block(QList<FieldInfo *> fieldInfos, QList<CouplingInfo*> couplings) :
    m_couplings(couplings)
{
    foreach (FieldInfo* fi, fieldInfos)
    {
        Field* field = new Field(fi);
        foreach (CouplingInfo* couplingInfo, Agros2D::problem()->couplingInfos())
        {
            if (couplingInfo->isWeak() && (couplingInfo->targetField() == fi))
            {
                field->m_couplingSources.push_back(couplingInfo);
            }
        }

        m_fields.append(field);
    }
}

Solver<double>* Block::prepareSolver()
{
    //Agros2D::log()->printDebug(QObject::tr("Solver"), QObject::tr("prepare solver"));

    Solver<double>* solver = new Solver<double>;

    foreach (Field* field, m_fields)
    {
        if (! field->solveInitVariables())
            assert(0); //TODO co to znamena?
    }

    solver->init(this);

    return solver;
}

bool Block::isTransient() const
{
    foreach (Field *field, m_fields)
    {
        if (field->fieldInfo()->analysisType() == AnalysisType_Transient)
            return true;
    }

    return false;
}

double Block::timeSkip() const
{
    double skip = 0.;
    foreach (Field *field, m_fields)
    {
        if (field->m_fieldInfo->analysisType() == AnalysisType_Transient)
            continue;

        double sActual = field->fieldInfo()->timeSkip().number();
        if ((skip == 0.) || (sActual < skip))
            skip = sActual;
    }
    return skip;
}

AdaptivityType Block::adaptivityType() const
{
    AdaptivityType at = m_fields.at(0)->fieldInfo()->adaptivityType();

    foreach (Field *field, m_fields)
    {
        // todo: ensure in GUI
        assert(field->fieldInfo()->adaptivityType() == at);
    }

    return at;
}

int Block::adaptivitySteps() const
{
    int as = m_fields.at(0)->fieldInfo()->adaptivitySteps();

    foreach (Field *field, m_fields)
    {
        // todo: ensure in GUI
        assert(field->fieldInfo()->adaptivitySteps() == as);
    }

    return as;
}

int Block::adaptivityBackSteps() const
{
    int abs = m_fields.at(0)->fieldInfo()->adaptivityBackSteps();

    foreach (Field *field, m_fields)
    {
        // todo: ensure in GUI
        assert(field->fieldInfo()->adaptivityBackSteps() == abs);
    }

    return abs;
}

int Block::adaptivityRedoneEach() const
{
    int re = m_fields.at(0)->fieldInfo()->adaptivityRedoneEach();

    foreach (Field *field, m_fields)
    {
        // todo: ensure in GUI
        assert(field->fieldInfo()->adaptivityRedoneEach() == re);
    }

    return re;
}

double Block::adaptivityTolerance() const
{
    double at = m_fields.at(0)->fieldInfo()->adaptivityTolerance();

    foreach (Field *field, m_fields)
    {
        // todo: ensure in GUI
        assert(field->fieldInfo()->adaptivityTolerance() == at);
    }

    return at;
}

int Block::numSolutions() const
{
    int num = 0;
    foreach (Field *field, m_fields)
        num += field->fieldInfo()->module()->numberOfSolutions();

    return num;
}

int Block::offset(Field *fieldParam) const
{
    int offset = 0;

    foreach (Field* field, m_fields)
    {
        if (field == fieldParam)
            return offset;
        else
            offset += field->fieldInfo()->module()->numberOfSolutions();
    }

    assert(0);
}

LinearityType Block::linearityType() const
{
    int linear = 0;
    int newton = 0;
    int picard = 0;

    foreach (Field* field, m_fields)
    {
        FieldInfo* fieldInfo = field->fieldInfo();
        if (fieldInfo->linearityType() == LinearityType_Linear)
            linear++;
        if (fieldInfo->linearityType() == LinearityType_Newton)
            newton++;
        if (fieldInfo->linearityType() == LinearityType_Picard)
            picard++;
    }
    assert(linear * newton * picard == 0); // all hard coupled fields has to be solved by the same method

    if (linear)
        return LinearityType_Linear;
    else if (newton)
        return  LinearityType_Newton;
    else if (picard)
        return  LinearityType_Picard;
    else
        assert(0);
}

double Block::nonlinearTolerance() const
{
    double tolerance = 10e20;

    foreach (Field* field, m_fields)
    {
        FieldInfo* fieldInfo = field->fieldInfo();
        if (fieldInfo->nonlinearTolerance() < tolerance)
            tolerance = fieldInfo->nonlinearTolerance();
    }

    return tolerance;
}

int Block::nonlinearSteps() const
{
    int steps = 0;

    foreach (Field* field, m_fields)
    {
        FieldInfo* fieldInfo = field->fieldInfo();
        if (fieldInfo->nonlinearSteps() > steps)
            steps = fieldInfo->nonlinearSteps();
    }

    return steps;
}

bool Block::newtonAutomaticDamping() const
{
    foreach (Field* field, m_fields)
    {
        FieldInfo* fieldInfo = field->fieldInfo();
        if (!fieldInfo->newtonAutomaticDamping())
            return false;
    }

    return true;
}

double Block::newtonDampingCoeff() const
{
    double coeff = 1.0;

    foreach (Field* field, m_fields)
    {
        FieldInfo* fieldInfo = field->fieldInfo();
        if (fieldInfo->newtonDampingCoeff() < coeff)
            coeff = fieldInfo->newtonDampingCoeff();
    }

    return coeff;
}

int Block::newtonDampingNumberToIncrease() const
{
    int number = 0;

    foreach (Field* field, m_fields)
    {
        FieldInfo* fieldInfo = field->fieldInfo();
        if (fieldInfo->newtonDampingNumberToIncrease() > number)
            number = fieldInfo->newtonDampingNumberToIncrease();
    }

    return number;
}

bool Block::picardAndersonAcceleration() const
{
    foreach (Field* field, m_fields)
    {
        FieldInfo* fieldInfo = field->fieldInfo();
        if (!fieldInfo->picardAndersonAcceleration())
            return false;
    }

    return true;
}

double Block::picardAndersonBeta() const
{
    double number = 0;

    foreach (Field* field, m_fields)
    {
        FieldInfo* fieldInfo = field->fieldInfo();
        // TODO: check ">"
        if (fieldInfo->picardAndersonBeta() > number)
            number = fieldInfo->picardAndersonBeta();
    }

    return number;
}

int Block::picardAndersonNumberOfLastVectors() const
{
    int number = 1;

    foreach (Field* field, m_fields)
    {
        FieldInfo* fieldInfo = field->fieldInfo();
        if (fieldInfo->picardAndersonNumberOfLastVectors() > number)
            number = fieldInfo->picardAndersonNumberOfLastVectors();
    }

    return number;
}

bool Block::contains(FieldInfo *fieldInfo) const
{
    foreach(Field* field, m_fields)
    {
        if (field->fieldInfo() == fieldInfo)
            return true;
    }
    return false;
}

Field* Block::field(FieldInfo *fieldInfo) const
{
    foreach (Field* field, m_fields)
    {
        if (fieldInfo == field->fieldInfo())
            return field;
    }

    return NULL;
}

Hermes::vector<Hermes::Hermes2D::ProjNormType> Block::projNormTypeVector() const
{
    Hermes::vector<Hermes::Hermes2D::ProjNormType> vec;

    foreach (Field* field, m_fields)
    {
        for (int comp = 0; comp < field->fieldInfo()->module()->numberOfSolutions(); comp++)
        {
            Hermes::Hermes2D::SpaceType spaceType = field->fieldInfo()->module()->spaces()[comp+1].type();
            if (spaceType == Hermes::Hermes2D::HERMES_H1_SPACE)
                vec.push_back(Hermes::Hermes2D::HERMES_H1_NORM);
            else if (spaceType == Hermes::Hermes2D::HERMES_L2_SPACE)
                vec.push_back(Hermes::Hermes2D::HERMES_L2_NORM);
            else
                assert(0);
        }
    }
    return vec;
}

ostream& operator<<(ostream& output, const Block& id)
{
    output << "Block ";
    return output;
}


