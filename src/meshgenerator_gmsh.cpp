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

#include "meshgenerator_gmsh.h"

#include "gui.h"
#include "scene.h"

#include "scenebasic.h"
#include "scenenode.h"
#include "sceneedge.h"
#include "scenelabel.h"

#include "sceneview_common.h"
#include "scenemarker.h"
#include "scenemarkerdialog.h"
#include "logview.h"

#include "hermes2d/module.h"
#include "hermes2d/module_agros.h"
#include "hermes2d/field.h"
#include "hermes2d/problem.h"

MeshGeneratorGMSH::MeshGeneratorGMSH() : MeshGenerator()
{    
}

bool MeshGeneratorGMSH::mesh()
{
    m_isError = false;

    QFile::remove(tempProblemFileName() + ".mesh");

    // create gmsh files
    if (writeToGmsh())
    {
        Util::log()->printDebug(tr("Mesh generator"), tr("Poly file was created"));

        // exec triangle
        QProcess processGmsh;
        processGmsh.setStandardOutputFile(tempProblemFileName() + ".gmsh.out");
        processGmsh.setStandardErrorFile(tempProblemFileName() + ".gmsh.err");
        connect(&processGmsh, SIGNAL(finished(int)), this, SLOT(meshGmshCreated(int)));

        QString gmshBinary = "gmsh";
        if (QFile::exists(QApplication::applicationDirPath() + QDir::separator() + "gmsh.exe"))
            gmshBinary = "\"" + QApplication::applicationDirPath() + QDir::separator() + "gmsh.exe\"";
        if (QFile::exists(QApplication::applicationDirPath() + QDir::separator() + "gmsh"))
            gmshBinary = QApplication::applicationDirPath() + QDir::separator() + "gmsh";

        processGmsh.start(QString(Util::config()->commandGmsh).
                          arg(gmshBinary).
                          arg(tempProblemFileName()));

        if (!processGmsh.waitForStarted(100000))
        {
            Util::log()->printError(tr("Mesh generator"), tr("could not start GMSH."));
            processGmsh.kill();

            m_isError = true;
        }
        else
        {
            // copy gmsh files
            if ((!Util::config()->deleteMeshFiles) && (!Util::problem()->config()->fileName().isEmpty()))
            {
                QFileInfo fileInfoOrig(Util::problem()->config()->fileName());

                QFile::copy(tempProblemFileName() + ".geo", fileInfoOrig.absolutePath() + QDir::separator() + fileInfoOrig.baseName() + ".geo");
            }

            while (!processGmsh.waitForFinished()) {}
        }
    }
    else
    {
        m_isError = true;
    }

    return !m_isError;
}

void MeshGeneratorGMSH::meshGmshCreated(int exitCode)
{
    if (exitCode == 0)
    {
        Util::log()->printMessage(tr("Mesh generator"), tr("mesh files were created"));
        // convert gmsh mesh to hermes mesh
        if (readGmshMeshFile())
        {
            Util::log()->printMessage(tr("Mesh generator"), tr("mesh was converted to Hermes2D mesh file"));

            // copy triangle files
            if ((!Util::config()->deleteHermes2DMeshFile) && (!Util::problem()->config()->fileName().isEmpty()))
            {
                QFileInfo fileInfoOrig(Util::problem()->config()->fileName());

                QFile::copy(tempProblemFileName() + ".mesh", fileInfoOrig.absolutePath() + QDir::separator() + fileInfoOrig.baseName() + ".mesh");
            }

            //  remove triangle temp files
            /*
            QFile::remove(tempProblemFileName() + ".geo");
            QFile::remove(tempProblemFileName() + ".msh");
            QFile::remove(tempProblemFileName() + ".gmsh.out");
            QFile::remove(tempProblemFileName() + ".gmsh.err");
            */
            Util::log()->printMessage(tr("Mesh generator"), tr("mesh files were deleted"));

            // load mesh
            try
            {
                QMap<FieldInfo*, Hermes::Hermes2D::Mesh*> meshes = readMeshesFromFile(tempProblemFileName() + ".xml");

                // FIXME: jinak
                Util::problem()->setMeshesInitial(meshes);
            }
            catch (Hermes::Exceptions::Exception& e)
            {
                m_isError = true;

                Util::log()->printError(tr("Mesh generator"), QString("%1").arg(e.what()));
            }

        }
        else
        {
            m_isError = true;
            // QFile::remove(Util::problem()->config()->fileName() + ".mesh");
        }
    }
    else
    {
        m_isError = true;
        QString errorMessage = readFileContent(Util::problem()->config()->fileName() + ".gmsh.out");
        Util::log()->printError(tr("Mesh generator"), errorMessage);
    }
}

bool MeshGeneratorGMSH::writeToGmsh()
{
    // basic check
    if (Util::scene()->nodes->length() < 3)
    {
        Util::log()->printError(tr("Mesh generator"), tr("invalid number of nodes (%1 < 3)").arg(Util::scene()->nodes->length()));
        return false;
    }
    if (Util::scene()->edges->length() < 3)
    {
        Util::log()->printError(tr("Mesh generator"), tr("invalid number of edges (%1 < 3)").arg(Util::scene()->edges->length()));
        return false;
    }

    // save current locale
    char *plocale = setlocale (LC_NUMERIC, "");
    setlocale (LC_NUMERIC, "C");

    QDir dir;
    dir.mkdir(QDir::temp().absolutePath() + "/agros2d");
    QFile file(tempProblemFileName() + ".geo");

    if (!file.open(QIODevice::WriteOnly | QIODevice::Text))
    {
        Util::log()->printError(tr("Mesh generator"), tr("could not create GMSH geo mesh file (%1)").arg(file.errorString()));
        return false;
    }
    QTextStream out(&file);

    out << QString("mesh_size = 0;\n");

    // nodes
    QString outNodes;
    int nodesCount = 0;
    for (int i = 0; i<Util::scene()->nodes->length(); i++)
    {
        outNodes += QString("Point(%1) = {%2, %3, 0, mesh_size};\n").
                arg(i).
                arg(Util::scene()->nodes->at(i)->point().x, 0, 'f', 10).
                arg(Util::scene()->nodes->at(i)->point().y, 0, 'f', 10);
        nodesCount++;
    }

    // edges
    QString outEdges;
    int edgesCount = 0;
    for (int i = 0; i<Util::scene()->edges->length(); i++)
    {
        if (Util::scene()->edges->at(i)->angle() == 0)
        {
            // line
            outEdges += QString("Line(%1) = {%2, %3};\n").
                    arg(edgesCount).
                    arg(Util::scene()->nodes->items().indexOf(Util::scene()->edges->at(i)->nodeStart())).
                    arg(Util::scene()->nodes->items().indexOf(Util::scene()->edges->at(i)->nodeEnd()));
            edgesCount++;
        }
        else
        {
            // arc
            // add pseudo nodes
            Point center = Util::scene()->edges->at(i)->center();
            outNodes += QString("Point(%1) = {%2, %3, 0};\n").
                    arg(nodesCount).
                    arg(center.x, 0, 'f', 10).
                    arg(center.y, 0, 'f', 10);
            nodesCount++;

            outEdges += QString("Circle(%1) = {%2, %3, %4};\n").
                    arg(edgesCount).
                    arg(Util::scene()->nodes->items().indexOf(Util::scene()->edges->at(i)->nodeStart())).
                    arg(nodesCount - 1).
                    arg(Util::scene()->nodes->items().indexOf(Util::scene()->edges->at(i)->nodeEnd()));

            edgesCount++;
        }
    }

    /*
    // holes
    int holesCount = 0;
    foreach (SceneLabel *label, Util::scene()->labels->items())
        if (label->markersCount() == 0)
            holesCount++;

    QString outHoles = QString("%1\n").arg(holesCount);
    holesCount = 0;
    foreach (SceneLabel *label, Util::scene()->labels->items())
    {
        if (label->markersCount() == 0)
        {
            outHoles += QString("%1  %2  %3\n").
                    arg(holesCount).
                    // arg(Util::scene()->labels->items().indexOf(label) + 1).
                    arg(label->point().x, 0, 'f', 10).
                    arg(label->point().y, 0, 'f', 10);

            holesCount++;
        }
    }

    // labels
    QString outLabels;
    int labelsCount = 0;
    foreach (SceneLabel *label, Util::scene()->labels->items())
    {
        if (label->markersCount() > 0)
        {
            outLabels += QString("%1  %2  %3  %4  %5\n").
                    arg(labelsCount).
                    arg(label->point().x, 0, 'f', 10).
                    arg(label->point().y, 0, 'f', 10).
                    // arg(labelsCount + 1). // triangle returns zero region number for areas without marker, markers must start from 1
                    arg(Util::scene()->labels->items().indexOf(label) + 1).
                    arg(label->area());
            labelsCount++;
        }
    }

    out << outHoles;
    outLabels.insert(0, QString("%1 1\n").
                     arg(labelsCount)); // - holes
    out << outLabels;
    */

    // TODO: find loops

    QString outLoops;
    outLoops.append(QString("Line Loop(1) = {0, 3, 2, 1};\n"));
    outLoops.append(QString("Plane Surface(1) = {1};\n"));
    outLoops.append(QString("Line Loop(2) = {4, 5, 6, -1};\n"));
    outLoops.append(QString("Plane Surface(2) = {2};\n"));
    outLoops.append("\n");

    // quad mesh
    if (Util::problem()->config()->meshType() == MeshType_GMSH_Quad)
        outLoops.append(QString("Recombine Surface {1, 2};\n"));

    // Mesh.Algorithm - 1=MeshAdapt, 2=Automatic, 5=Delaunay, 6=Frontal, 7=bamg, 8=delquad
    QString outCommands;
    if (Util::problem()->config()->meshType() == MeshType_GMSH_Quad)
    {
        outCommands.append(QString("Mesh.Algorithm = 8;\n"));
        outCommands.append(QString("Mesh.SubdivisionAlgorithm = 1;\n"));
    }
    else if (Util::problem()->config()->meshType() == MeshType_GMSH_Triangle)
    {
        outCommands.append(QString("Mesh.Algorithm = 2;\n"));
    }

    outNodes.insert(0, QString("\n// nodes\n"));
    out << outNodes;
    outEdges.insert(0, QString("\n// edges\n"));
    out << outEdges;
    outLoops.insert(0, QString("\n// loops\n"));
    out << outLoops;
    outCommands.insert(0, QString("\n// commands\n"));
    out << outCommands;

    file.waitForBytesWritten(0);
    file.close();

    // set system locale
    setlocale(LC_NUMERIC, plocale);

    return true;
}

bool MeshGeneratorGMSH::readGmshMeshFile()
{
    nodeList.clear();
    edgeList.clear();
    elementList.clear();

    int k;

    QFile fileGMSH(tempProblemFileName() + ".msh");
    if (!fileGMSH.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        Util::log()->printError(tr("Mesh generator"), tr("could not read GMSH mesh file"));
        return false;
    }
    QTextStream inGMSH(&fileGMSH);

    // nodes
    inGMSH.readLine();
    inGMSH.readLine();
    inGMSH.readLine();
    inGMSH.readLine();
    sscanf(inGMSH.readLine().toStdString().c_str(), "%i", &k);
    for (int i = 0; i < k; i++)
    {
        int n;
        double x, y, z;

        sscanf(inGMSH.readLine().toStdString().c_str(), "%i %lf %lf %lf", &n, &x, &y, &z);
        nodeList.append(Point(x, y));
    }

    // elements
    inGMSH.readLine();
    inGMSH.readLine();
    sscanf(inGMSH.readLine().toStdString().c_str(), "%i", &k);
    QSet<int> labelMarkersCheck;
    for (int i = 0; i < k; i++)
    {
        int quad[4];
        int n, type, phys, part, marker;

        if (sscanf(inGMSH.readLine().toStdString().c_str(), "%i %i %i %i %i %i %i %i %i",
                   &n, &type, &phys, &part, &marker, &quad[0], &quad[1], &quad[2], &quad[3]))
        {
            // edge
            if (type == 1)
                edgeList.append(MeshEdge(quad[0] - 1, quad[1] - 1, marker)); // marker conversion from gmsh, where it starts from 1
            // triangle
            if (type == 2)
                elementList.append(MeshElement(quad[0] - 1, quad[1] - 1, quad[2] - 1, marker - 1)); // marker conversion from gmsh, where it starts from 1
            // quad
            if (type == 3)
                elementList.append(MeshElement(quad[0] - 1, quad[1] - 1, quad[2] - 1, quad[3] - 1, marker - 1)); // marker conversion from gmsh, where it starts from 1
        }
        /*

        if (marker == 0)
        {
            Util::log()->printError(tr("Mesh generator"), tr("some areas have no label marker"));
            return false;
        }
        */
        labelMarkersCheck.insert(marker - 1);
    }
    int elementCountLinear = elementList.count();

    fileGMSH.close();

    writeToHermes();

    nodeList.clear();
    edgeList.clear();
    elementList.clear();  
}
