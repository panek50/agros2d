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

#include "sceneview_data.h"
#include "sceneview.h"

// scene view
static SceneView *m_sceneView = NULL;;

SceneView *sceneView()
{
    return m_sceneView;
}

static inline double* computeNormal(double p0x, double p0y, double p0z, double p1x, double p1y, double p1z, double p2x, double p2y, double p2z)
{
    double ax = (p1x - p0x);
    double ay = (p1y - p0y);
    double az = (p1z - p0z);

    double bx = (p2x - p0x);
    double by = (p2y - p0y);
    double bz = (p2z - p0z);

    double nx = ay * bz - az * by;
    double ny = az * bx - ax * bz;
    double nz = ax * by - ay * bx;

    // normalize
    // double l = 1.0 / sqrt(sqr(nx) + sqr(ny) + sqr(nz));
    // double p[3] = { nx*l, ny*l, nz*l };

    double p[3] = { nx, ny, nz };

    return p;
}

// *******************************************************************************************************

SceneViewSettings::SceneViewSettings()
{
    defaultValues();
}

void SceneViewSettings::defaultValues()
{
    scalarRangeMin =  CONST_DOUBLE;
    scalarRangeMax = -CONST_DOUBLE;

    // visible objects
    showGeometry = true;
    showGrid = true;
    showInitialMesh = false;

    postprocessorShow = SceneViewPostprocessorShow_ScalarView;

    showContours = false;
    showVectors = false;
    showSolutionMesh = false;

    contourPhysicFieldVariable = Util::scene()->problemInfo()->hermes()->contourPhysicFieldVariable();

    scalarPhysicFieldVariable = Util::scene()->problemInfo()->hermes()->scalarPhysicFieldVariable();
    scalarPhysicFieldVariableComp = Util::scene()->problemInfo()->hermes()->scalarPhysicFieldVariableComp();
    scalarRangeAuto = true;

    vectorPhysicFieldVariable = Util::scene()->problemInfo()->hermes()->vectorPhysicFieldVariable();
}

// *******************************************************************************************************

SceneView::SceneView(QWidget *parent): QGLWidget(QGLFormat(QGL::SampleBuffers), parent)
{
    m_sceneView = this;
    m_mainWindow = (QMainWindow *) parent;
    m_scene = Util::scene();

    connect(m_scene->sceneSolution(), SIGNAL(timeStepChanged(bool)), this, SLOT(timeStepChanged(bool)));
    connect(m_scene->sceneSolution(), SIGNAL(solved()), this, SLOT(solved()));
    connect(m_scene->sceneSolution(), SIGNAL(processedRangeContour()), this, SLOT(processedRangeContour()));
    connect(m_scene->sceneSolution(), SIGNAL(processedRangeScalar()), this, SLOT(processedRangeScalar()));
    connect(m_scene->sceneSolution(), SIGNAL(processedRangeVector()), this, SLOT(processedRangeVector()));
    connect(m_scene->sceneSolution(), SIGNAL(meshed()), this, SLOT(clearGLLists()));

    connect(m_scene, SIGNAL(invalidated()), this, SLOT(doInvalidated()));
    connect(m_scene, SIGNAL(defaultValues()), this, SLOT(doDefaultValues()));

    createActions();
    createMenu();

    doDefaultValues();

    setMouseTracking(true);
    setFocusPolicy(Qt::StrongFocus);
    setContextMenuPolicy(Qt::DefaultContextMenu);

#ifdef Q_WS_X11
    setFont(QFont("Monospace", 9));
#endif
#ifdef Q_WS_WIN
    setFont(QFont("Courier New", 9));
#endif
#ifdef Q_WS_MAC
    setFont(QFont("Monaco", 12));
#endif

    setMinimumSize(400, 400);
}

SceneView::~SceneView()
{
}

void SceneView::createActions()
{
    // scene - zoom
    actSceneZoomIn = new QAction(icon("zoom-in"), tr("Zoom in"), this);
    actSceneZoomIn->setShortcut(QKeySequence::ZoomIn);
    actSceneZoomIn->setStatusTip(tr("Zoom in"));
    connect(actSceneZoomIn, SIGNAL(triggered()), this, SLOT(doZoomIn()));

    actSceneZoomOut = new QAction(icon("zoom-out"), tr("Zoom out"), this);
    actSceneZoomOut->setShortcut(QKeySequence::ZoomOut);
    actSceneZoomOut->setStatusTip(tr("Zoom out"));
    connect(actSceneZoomOut, SIGNAL(triggered()), this, SLOT(doZoomOut()));

    actSceneZoomBestFit = new QAction(icon("zoom-original"), tr("Zoom best fit"), this);
    actSceneZoomBestFit->setStatusTip(tr("Best fit"));
    connect(actSceneZoomBestFit, SIGNAL(triggered()), this, SLOT(doZoomBestFit()));

    actSceneZoomRegion = new QAction(icon("zoom-best-fit"), tr("Zoom region"), this);
    actSceneZoomRegion->setStatusTip(tr("Zoom region"));
    actSceneZoomRegion->setCheckable(true);

    // scene - operate on items
    actSceneModeNode = new QAction(icon("scene-node"), tr("Operate on &nodes"), this);
    actSceneModeNode->setShortcut(Qt::Key_F5);
    actSceneModeNode->setStatusTip(tr("Operate on nodes"));
    actSceneModeNode->setCheckable(true);

    actSceneModeEdge = new QAction(icon("scene-edge"), tr("Operate on &edges"), this);
    actSceneModeEdge->setShortcut(Qt::Key_F6);
    actSceneModeEdge->setStatusTip(tr("Operate on edges"));
    actSceneModeEdge->setCheckable(true);

    actSceneModeLabel = new QAction(icon("scene-label"), tr("Operate on &labels"), this);
    actSceneModeLabel->setShortcut(Qt::Key_F7);
    actSceneModeLabel->setStatusTip(tr("Operate on labels"));
    actSceneModeLabel->setCheckable(true);

    actSceneModePostprocessor = new QAction(icon("scene-postprocessor"), tr("&Postprocessor"), this);
    actSceneModePostprocessor->setShortcut(Qt::Key_F8);
    actSceneModePostprocessor->setStatusTip(tr("Postprocessor"));
    actSceneModePostprocessor->setCheckable(true);

    actSceneModeGroup = new QActionGroup(this);
    actSceneModeGroup->addAction(actSceneModeNode);
    actSceneModeGroup->addAction(actSceneModeEdge);
    actSceneModeGroup->addAction(actSceneModeLabel);
    actSceneModeGroup->addAction(actSceneModePostprocessor);
    connect(actSceneModeGroup, SIGNAL(triggered(QAction *)), this, SLOT(doSceneModeSet(QAction *)));

    // material
    actMaterialGroup = new QActionGroup(this);
    connect(actMaterialGroup, SIGNAL(triggered(QAction *)), this, SLOT(doMaterialGroup(QAction *)));

    // boundary
    actBoundaryGroup = new QActionGroup(this);
    connect(actBoundaryGroup, SIGNAL(triggered(QAction *)), this, SLOT(doBoundaryGroup(QAction *)));

    // show
    actShowSolutionMesh = new QAction(tr("Solution mesh"), this);
    actShowSolutionMesh->setCheckable(true);

    actShowContours = new QAction(tr("Contours"), this);
    actShowContours->setCheckable(true);

    actShowVectors = new QAction(tr("Vectors"), this);
    actShowVectors->setCheckable(true);

    actShowGroup = new QActionGroup(this);
    actShowGroup->setExclusive(false);
    connect(actShowGroup, SIGNAL(triggered(QAction *)), this, SLOT(doShowGroup(QAction *)));
    actShowGroup->addAction(actShowSolutionMesh);
    actShowGroup->addAction(actShowContours);
    actShowGroup->addAction(actShowVectors);

    // postprocessor group
    actPostprocessorModeLocalPointValue = new QAction(icon("mode-localpointvalue"), tr("Local Values"), this);
    actPostprocessorModeLocalPointValue->setCheckable(true);

    actPostprocessorModeSurfaceIntegral = new QAction(icon("mode-surfaceintegral"), tr("Surface Integrals"), this);
    actPostprocessorModeSurfaceIntegral->setCheckable(true);

    actPostprocessorModeVolumeIntegral = new QAction(icon("mode-volumeintegral"), tr("Volume Integrals"), this);
    actPostprocessorModeVolumeIntegral->setCheckable(true);

    actPostprocessorModeGroup = new QActionGroup(this);
    connect(actPostprocessorModeGroup, SIGNAL(triggered(QAction *)), this, SLOT(doPostprocessorModeGroup(QAction*)));
    actPostprocessorModeGroup->addAction(actPostprocessorModeLocalPointValue);
    actPostprocessorModeGroup->addAction(actPostprocessorModeSurfaceIntegral);
    actPostprocessorModeGroup->addAction(actPostprocessorModeVolumeIntegral);

    // scene properties
    actSceneViewProperties = new QAction(icon("scene-properties"), tr("&Scene properties"), this);
    actSceneViewProperties->setShortcut(Qt::Key_F12);
    connect(actSceneViewProperties, SIGNAL(triggered()), this, SLOT(doSceneViewProperties()));

    // object properties
    actSceneObjectProperties = new QAction(icon("scene-properties"), tr("Object properties"), this);
    connect(actSceneObjectProperties, SIGNAL(triggered()), this, SLOT(doSceneObjectProperties()));

    // select region
    actSceneViewSelectRegion = new QAction(icon("scene-select-region"), tr("&Select region"), this);
    actSceneViewSelectRegion->setStatusTip(tr("Select region"));
    actSceneViewSelectRegion->setCheckable(true);

    actSceneViewSelectMarker = new QAction(icon(""), tr("Select by marker"), this);
    actSceneViewSelectMarker->setStatusTip(tr("Select by marker"));
    connect(actSceneViewSelectMarker, SIGNAL(triggered()), this, SLOT(doSelectMarker()));
}

void SceneView::createMenu()
{
    mnuInfo = new QMenu(this);

    QMenu *mnuModeGroup = new QMenu(tr("Set mode"), this);
    mnuModeGroup->addAction(actSceneModeNode);
    mnuModeGroup->addAction(actSceneModeEdge);
    mnuModeGroup->addAction(actSceneModeLabel);
    mnuModeGroup->addAction(actSceneModePostprocessor);

    mnuInfo->addAction(m_scene->actNewNode);
    mnuInfo->addAction(m_scene->actNewEdge);
    mnuInfo->addAction(m_scene->actNewLabel);
    mnuInfo->addSeparator();
    mnuInfo->addAction(m_scene->actNewEdgeMarker);
    mnuInfo->addAction(m_scene->actNewLabelMarker);
    mnuInfo->addSeparator();
    mnuInfo->addAction(actSceneViewSelectRegion);
    mnuInfo->addAction(m_scene->actTransform);
    mnuInfo->addSeparator();
    mnuInfo->addMenu(mnuModeGroup);
    mnuInfo->addSeparator();
    mnuInfo->addAction(actSceneObjectProperties);
    mnuInfo->addAction(m_scene->actProblemProperties);
    mnuInfo->addAction(actSceneViewProperties);
}

void SceneView::initializeGL()
{
    glShadeModel(GL_SMOOTH);    
    glEnable(GL_NORMALIZE);

    clearGLLists();
}

void SceneView::resizeGL(int w, int h)
{
    setupViewport(w, h);

    if (Util::scene()->sceneSolution()->isSolved() && m_sceneMode == SceneMode_Postprocessor)
    {
        paletteFilter();
        paletteUpdateTexAdjust();
        paletteCreate();
    }
}

void SceneView::loadProjection2d(bool setScene)
{
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(3.0, contextWidth()-6.0, contextHeight()-6.0, 3.0, -10.0, -10.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (setScene)
    {
        glScaled(m_scale2d/aspect(), m_scale2d, m_scale2d);

        glTranslated(-m_offset2d.x, -m_offset2d.y, 0.0);
    }
}

void SceneView::loadProjection3d(bool setScene)
{
    int fov = 50.0;
    double znear = 0.001;
    double zfar = 100.0;

    double right = znear * tan((double) fov / 2.0 / 180.0 * M_PI);
    double top = (double) contextHeight() / contextWidth() * right;
    double left = -right;
    double bottom = -top;
    double offsx = (right - left) / contextWidth();
    double offsy = (top - bottom) / contextHeight();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    // glFrustum(left - offsx, right - offsx, bottom - offsy, top - offsy, znear, zfar);
    glOrtho(3.0, contextWidth()-6.0, contextHeight()-6.0, 3.0, -10.0, -10.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    if (setScene)
    {        
        glScaled(1.0/aspect(), 1.0, 1.0);

        // move to origin
        RectPoint rect = Util::scene()->boundingBox();
        glTranslated(-m_offset3d.x, -m_offset3d.y, 0.0);

        glRotated(m_rotation3d.x, 1.0, 0.0, 0.0);
        glRotated(m_rotation3d.z, 0.0, 1.0, 0.0);
        glRotated(m_rotation3d.y, 0.0, 0.0, 1.0);

        if (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_ScalarView3D)
        {
            glTranslated(- m_scale3d * (rect.start.x + rect.end.x) / 2.0, - m_scale3d * (rect.start.y + rect.end.y) / 2.0, 0.0);
        }
        else
        {
            if (Util::scene()->problemInfo()->problemType == ProblemType_Planar)
            {
                glTranslated(- m_scale3d * (rect.start.x + rect.end.x) / 2.0, - m_scale3d * (rect.start.y + rect.end.y) / 2.0, 0.0);
            }
            else
            {
                glTranslated(0.0, - m_scale3d * (rect.start.y + rect.end.y) / 2.0, 0.0);
            }
        }

        glScaled(m_scale3d, m_scale3d, m_scale3d);
    }
}

void SceneView::setupViewport(int w, int h)
{
    glViewport(0, 0, w, h);
}

QPixmap SceneView::renderScenePixmap(int w, int h, bool useContext)
{
    QPixmap pixmap = renderPixmap(w, h, useContext);

    resizeGL(contextWidth(), contextHeight());

    return pixmap;
}

void SceneView::paintGL()
{
    glClearColor(Util::config()->colorBackground.redF(),
                 Util::config()->colorBackground.greenF(),
                 Util::config()->colorBackground.blueF(), 0);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // if (Util::scene()->sceneSolution()->isSolving())
    //    return;

    if (is3DMode())
    {
        glClear(GL_DEPTH_BUFFER_BIT);

        if (m_scene->sceneSolution()->isMeshed() && (m_sceneMode == SceneMode_Postprocessor))
        {
            if (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_Model) paintScalarField3DSolid();
        }

        if (m_scene->sceneSolution()->isSolved() && (m_sceneMode == SceneMode_Postprocessor))
        {
            if (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_ScalarView3D) paintScalarField3D();
            if (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_ScalarView3DSolid) paintScalarField3DSolid();

            if (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_ScalarView3D ||
                m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_ScalarView3DSolid)
                paintScalarFieldColorBar(m_sceneViewSettings.scalarRangeMin, m_sceneViewSettings.scalarRangeMax);
        }
    }
    else
    {
        glDisable(GL_DEPTH_TEST);

        // grid
        if (m_sceneViewSettings.showGrid) paintGrid();

        // view
        if (m_scene->sceneSolution()->isSolved() && (m_sceneMode == SceneMode_Postprocessor))
        {
            if (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_ScalarView) paintScalarField();
            if (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_Order) paintOrder();

            if (m_sceneViewSettings.showContours) paintContours();
            if (m_sceneViewSettings.showVectors) paintVectors();
            if (m_sceneViewSettings.showSolutionMesh) paintSolutionMesh();
        }

        // initial mesh
        if (m_sceneViewSettings.showInitialMesh) paintInitialMesh();

        // geometry
        if (m_sceneViewSettings.showGeometry) paintGeometry();

        if (m_scene->sceneSolution()->isSolved() && (m_sceneMode == SceneMode_Postprocessor))
        {
            if (actPostprocessorModeVolumeIntegral->isChecked()) paintPostprocessorSelectedVolume();
            if (actPostprocessorModeSurfaceIntegral->isChecked()) paintPostprocessorSelectedSurface();

            // bars
            if (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_ScalarView) paintScalarFieldColorBar(m_sceneViewSettings.scalarRangeMin, m_sceneViewSettings.scalarRangeMax);
            if (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_Order) paintOrderColorBar();
        }

        // rulers
        if (Util::config()->showRulers) paintRulers();

        paintZoomRegion();
        paintSnapToGrid();
        paintChartLine();
    }

    paintSceneModeLabel();
}

void SceneView::clearGLLists()
{
    if (m_listContours != -1) glDeleteLists(m_listContours, 1);
    if (m_listVectors != -1) glDeleteLists(m_listVectors, 1);
    if (m_listScalarField != -1) glDeleteLists(m_listScalarField, 1);
    if (m_listScalarField3D != -1) glDeleteLists(m_listScalarField3D, 1);
    if (m_listScalarField3DSolid != -1) glDeleteLists(m_listScalarField3DSolid, 1);
    if (m_listOrder != -1) glDeleteLists(m_listOrder, 1);
    if (m_listModel != -1) glDeleteLists(m_listModel, 1);

    m_listContours = -1;
    m_listVectors = -1;
    m_listScalarField = -1;
    m_listScalarField3D = -1;
    m_listScalarField3DSolid = -1;
    m_listOrder = -1;
    m_listModel = -1;
}

// paint *****************************************************************************************************************************

void SceneView::paintBackground()
{
    // background
    glPushMatrix();
    glLoadIdentity();

    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();

    glOrtho(-1.0, 1.0, -1.0, 1.0, -10.0, -10.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glBegin(GL_QUADS);
    if (Util::config()->scalarView3DBackground)
        glColor3f(0.99, 0.99, 0.99);
    else
        glColor3f(Util::config()->colorBackground.redF(),
                  Util::config()->colorBackground.greenF(),
                  Util::config()->colorBackground.blueF());
    glVertex3d(-1.0, -1.0, 0.0);
    glVertex3d(1.0, -1.0, 0.0);
    if (Util::config()->scalarView3DBackground)
        glColor3f(0.44, 0.56, 0.89);
    glVertex3d(1.0, 1.0, 0.0);
    glVertex3d(-1.0, 1.0, 0.0);
    glEnd();

    glDisable(GL_POLYGON_OFFSET_FILL);

    glPopMatrix();
}

void SceneView::paintGrid()
{
    loadProjection2d(true);

    Point cornerMin = position(Point(0, 0));
    Point cornerMax = position(Point(contextWidth(), contextHeight()));

    glDisable(GL_DEPTH_TEST);

    int step = (((int) ((cornerMax - cornerMin).x / Util::config()->gridStep) + 1) / 5);
    if (step > 0.0)
    {
        glColor3f(Util::config()->colorGrid.redF(),
                  Util::config()->colorGrid.greenF(),
                  Util::config()->colorGrid.blueF());
        glLineWidth(1.0);
        glEnable(GL_LINE_STIPPLE);
        glLineStipple(1, 0x1C47);
        glBegin(GL_LINES);

        if ((((cornerMax.x-cornerMin.x)/Util::config()->gridStep + (cornerMin.y-cornerMax.y)/Util::config()->gridStep) < 200) &&
             ((cornerMax.x-cornerMin.x)/Util::config()->gridStep > 0) && ((cornerMin.y-cornerMax.y)/Util::config()->gridStep > 0))
        {
            // vertical lines
            for (int i = 0; i<cornerMax.x/Util::config()->gridStep; i++)
            {
                if ((step > 0) && i % step == 0)
                    glColor3f(Util::config()->colorCross.redF(),
                              Util::config()->colorCross.greenF(),
                              Util::config()->colorCross.blueF());
                else
                    glColor3f(Util::config()->colorGrid.redF(),
                              Util::config()->colorGrid.greenF(),
                              Util::config()->colorGrid.blueF());
                glVertex2d(i*Util::config()->gridStep, cornerMin.y);
                glVertex2d(i*Util::config()->gridStep, cornerMax.y);
            }
            for (int i = 0; i>cornerMin.x/Util::config()->gridStep; i--)
            {
                if ((step > 0) && i % step == 0)
                    glColor3f(Util::config()->colorCross.redF(),
                              Util::config()->colorCross.greenF(),
                              Util::config()->colorCross.blueF());
                else
                    glColor3f(Util::config()->colorGrid.redF(),
                              Util::config()->colorGrid.greenF(),
                              Util::config()->colorGrid.blueF());
                glVertex2d(i*Util::config()->gridStep, cornerMin.y);
                glVertex2d(i*Util::config()->gridStep, cornerMax.y);
            }

            // horizontal lines
            for (int i = 0; i<cornerMin.y/Util::config()->gridStep; i++)
            {
                if ((step > 0) && i % step == 0)
                    glColor3f(Util::config()->colorCross.redF(),
                              Util::config()->colorCross.greenF(),
                              Util::config()->colorCross.blueF());
                else
                    glColor3f(Util::config()->colorGrid.redF(),
                              Util::config()->colorGrid.greenF(),
                              Util::config()->colorGrid.blueF());
                glVertex2d(cornerMin.x, i*Util::config()->gridStep);
                glVertex2d(cornerMax.x, i*Util::config()->gridStep);
            }
            for (int i = 0; i>cornerMax.y/Util::config()->gridStep; i--)
            {
                if ((step > 0) && i % step == 0)
                    glColor3f(Util::config()->colorCross.redF(),
                              Util::config()->colorCross.greenF(),
                              Util::config()->colorCross.blueF());
                else
                    glColor3f(Util::config()->colorGrid.redF(),
                              Util::config()->colorGrid.greenF(),
                              Util::config()->colorGrid.blueF());
                glVertex2d(cornerMin.x, i*Util::config()->gridStep);
                glVertex2d(cornerMax.x, i*Util::config()->gridStep);
            }
        }
        glEnd();
        glDisable(GL_LINE_STIPPLE);
    }

    // axes
    glColor3f(Util::config()->colorCross.redF(),
              Util::config()->colorCross.greenF(),
              Util::config()->colorCross.blueF());
    glLineWidth(1.0);
    glBegin(GL_LINES);
    // y axis
    glVertex2d(0, cornerMin.y);
    glVertex2d(0, cornerMax.y);
    // x axis
    glVertex2d(cornerMin.x, 0);
    glVertex2d(cornerMax.x, 0);
    glEnd();
}

void SceneView::paintRulers()
{
    loadProjection2d(true);

    Point cornerMin = position(Point(0, 0));
    Point cornerMax = position(Point(contextWidth(), contextHeight()));

    // rulers
    double step = (((int) ((cornerMax - cornerMin).x / Util::config()->gridStep) + 1) / 5) * Util::config()->gridStep;

    Point size((2.0/contextWidth()*fontMetrics().width(" "))/m_scale2d*aspect(),
               (2.0/contextHeight()*fontMetrics().height())/m_scale2d);

    if (step > 0.0)
    {
        QString text;

        if (((cornerMax.x-cornerMin.x)/step > 0) && ((cornerMin.y-cornerMax.y)/step > 0))
        {
            glColor3f(0.3, 0.2, 0.0);

            // horizontal ticks
            for (int i = 0; i<cornerMax.x/step; i++)
            {
                text = QString::number(i*step, 'g', 4);
                renderTextPos(i*step - size.x*text.size() / 2.0, cornerMax.y + size.x / 4.0, 0.0, text);
            }
            for (int i = 0; i>cornerMin.x/step; i--)
            {
                text = QString::number(i*step, 'g', 4);
                renderTextPos(i*step - size.x*text.size() / 2.0, cornerMax.y + size.x / 4.0, 0.0, text);
            }

            // vertical ticks
            for (int i = 0; i<cornerMin.y/step; i++)
            {
                text = QString::number(i*step, 'g', 4);
                renderTextPos(cornerMin.x + size.x/4.0, i*step - size.y / 4.0, 0.0, text);
            }
            for (int i = 0; i>cornerMax.y/step; i--)
            {
                text = QString::number(i*step, 'g', 4);
                renderTextPos(cornerMin.x + size.x/4.0, i*step - size.y / 4.0, 0.0, text);
            }
        }
    }
}

void SceneView::paintGeometry()
{
    loadProjection2d(true);

    // edges
    foreach (SceneEdge *edge, m_scene->edges)
    {
        // edge with marker
        if (m_sceneMode == SceneMode_OperateOnEdges && edge->marker->type == PhysicFieldBC_None)
        {
            glEnable(GL_LINE_STIPPLE);
            glLineStipple(1, 0x8FFF);
        }

        glColor3f(Util::config()->colorEdges.redF(),
                  Util::config()->colorEdges.greenF(),
                  Util::config()->colorEdges.blueF());
        glLineWidth(Util::config()->edgeWidth);
        if (edge->isHighlighted)
        {
            glColor3f(Util::config()->colorHighlighted.redF(),
                      Util::config()->colorHighlighted.greenF(),
                      Util::config()->colorHighlighted.blueF());
            glLineWidth(Util::config()->edgeWidth + 2.0);
        }
        if (edge->isSelected)
        {
            glColor3f(Util::config()->colorSelected.redF(),
                      Util::config()->colorSelected.greenF(),
                      Util::config()->colorSelected.blueF());
            glLineWidth(Util::config()->edgeWidth + 2.0);
        }

        if (edge->angle == 0)
        {
            glBegin(GL_LINES);
            glVertex2d(edge->nodeStart->point.x, edge->nodeStart->point.y);
            glVertex2d(edge->nodeEnd->point.x, edge->nodeEnd->point.y);
            glEnd();
        }
        else
        {
            Point center = edge->center();
            double radius = edge->radius();
            double startAngle = atan2(center.y - edge->nodeStart->point.y, center.x - edge->nodeStart->point.x) / M_PI*180.0 - 180.0;

            drawArc(center, radius, startAngle, edge->angle, edge->angle/2.0);
        }

        glDisable(GL_LINE_STIPPLE);
        glLineWidth(1.0);
    }

    // nodes
    if (!(m_sceneMode == SceneMode_Postprocessor))
    {
        foreach (SceneNode *node, m_scene->nodes)
        {
            glColor3f(Util::config()->colorNodes.redF(),
                      Util::config()->colorNodes.greenF(),
                      Util::config()->colorNodes.blueF());
            glPointSize(Util::config()->nodeSize);

            glBegin(GL_POINTS);
            glVertex2d(node->point.x, node->point.y);
            glEnd();

            glColor3f(Util::config()->colorBackground.redF(),
                      Util::config()->colorBackground.greenF(),
                      Util::config()->colorBackground.blueF());
            glPointSize(Util::config()->nodeSize - 2.0);

            glBegin(GL_POINTS);
            glVertex2d(node->point.x, node->point.y);
            glEnd();

            if ((node->isSelected) || (node->isHighlighted))
            {
                if (node->isHighlighted)
                    glColor3f(Util::config()->colorHighlighted.redF(),
                              Util::config()->colorHighlighted.greenF(),
                              Util::config()->colorHighlighted.blueF());
                if (node->isSelected)
                    glColor3f(Util::config()->colorSelected.redF(),
                              Util::config()->colorSelected.greenF(),
                              Util::config()->colorSelected.blueF());

                glPointSize(Util::config()->nodeSize - 2.0);
                glBegin(GL_POINTS);
                glVertex2d(node->point.x, node->point.y);
                glEnd();
            }
        }

        glLineWidth(1.0);
    }
    // labels
    if (!(m_sceneMode == SceneMode_Postprocessor))
    {
        foreach (SceneLabel *label, m_scene->labels)
        {
            glColor3f(Util::config()->colorLabels.redF(),
                      Util::config()->colorLabels.greenF(),
                      Util::config()->colorLabels.blueF());
            glPointSize(Util::config()->labelSize);
            glBegin(GL_POINTS);
            glVertex2d(label->point.x, label->point.y);
            glEnd();

            glColor3f(Util::config()->colorBackground.redF(),
                      Util::config()->colorBackground.greenF(),
                      Util::config()->colorBackground.blueF());
            glPointSize(Util::config()->labelSize - 2.0);
            glBegin(GL_POINTS);
            glVertex2d(label->point.x, label->point.y);
            glEnd();

            if ((label->isSelected) || (label->isHighlighted))
            {
                if (label->isHighlighted)
                    glColor3f(Util::config()->colorHighlighted.redF(),
                              Util::config()->colorHighlighted.greenF(),
                              Util::config()->colorHighlighted.blueF());
                if (label->isSelected)
                    glColor3f(Util::config()->colorSelected.redF(),
                              Util::config()->colorSelected.greenF(),
                              Util::config()->colorSelected.blueF());

                glPointSize(Util::config()->labelSize - 2.0);
                glBegin(GL_POINTS);
                glVertex2d(label->point.x, label->point.y);
                glEnd();
            }
            glLineWidth(1.0);

            if (m_sceneMode == SceneMode_OperateOnLabels)
            {
                glColor3f(0.1, 0.1, 0.1);

                Point point;
                point.x = 2.0/contextWidth()*aspect()*fontMetrics().width(label->marker->name)/m_scale2d/2.0;
                point.y = 2.0/contextHeight()*fontMetrics().height()/m_scale2d;

                renderTextPos(label->point.x-point.x, label->point.y-point.y, 0.0, label->marker->name, false);
            }

            // area size
            if ((m_sceneMode == SceneMode_OperateOnLabels) || (m_sceneViewSettings.showInitialMesh))
            {
                double radius = sqrt(label->area/M_PI);
                glColor3f(0, 0.95, 0.9);
                glBegin(GL_LINE_LOOP);
                for (int i = 0; i<360; i = i + 10)
                {
                    glVertex2d(label->point.x + radius*cos(i/180.0*M_PI), label->point.y + radius*sin(i/180.0*M_PI));
                }
                glEnd();
            }
        }
    }
}

void SceneView::paintInitialMesh()
{
    if (!m_scene->sceneSolution()->isMeshed()) return;

    loadProjection2d(true);

    // draw initial mesh    
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(Util::config()->colorInitialMesh.redF(),
              Util::config()->colorInitialMesh.greenF(),
              Util::config()->colorInitialMesh.blueF());
    glLineWidth(1.3);

    // triangles
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < m_scene->sceneSolution()->meshInitial()->get_num_active_elements(); i++)
    {
        Element *element = m_scene->sceneSolution()->meshInitial()->get_element_fast(i);
        if (element->is_triangle())
        {
            glVertex2d(element->vn[0]->x, element->vn[0]->y);
            glVertex2d(element->vn[1]->x, element->vn[1]->y);
            glVertex2d(element->vn[2]->x, element->vn[2]->y);
        }
    }
    glEnd();
}

void SceneView::paintSolutionMesh()
{
    if (!m_isSolutionPrepared) return;

    loadProjection2d(true);

    // draw solution mesh
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
    glColor3f(Util::config()->colorSolutionMesh.redF(),
              Util::config()->colorSolutionMesh.greenF(),
              Util::config()->colorSolutionMesh.blueF());
    glLineWidth(1.0);

    // triangles
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < m_scene->sceneSolution()->sln()->get_mesh()->get_num_elements(); i++)
    {
        Element *element = m_scene->sceneSolution()->sln()->get_mesh()->get_element_fast(i);
        if (element->is_triangle())
        {
            glVertex2d(element->vn[0]->x, element->vn[0]->y);
            glVertex2d(element->vn[1]->x, element->vn[1]->y);
            glVertex2d(element->vn[2]->x, element->vn[2]->y);
        }
    }
    glEnd();
}

void SceneView::paintOrder()
{
    if (!m_isSolutionPrepared) return;

    loadProjection2d(true);

    if (m_listOrder == -1)
    {
        m_listOrder = glGenLists(1);
        glNewList(m_listOrder, GL_COMPILE);

        // order scalar view
        m_scene->sceneSolution()->ordView().lock_data();

        double3* vert = m_scene->sceneSolution()->ordView().get_vertices();
        int3* tris = m_scene->sceneSolution()->ordView().get_triangles();

        // draw mesh
        int min = 11;
        int max = 1;
        for (int i = 0; i < m_scene->sceneSolution()->ordView().get_num_triangles(); i++)
        {
            if (vert[tris[i][0]][2] < min) min = vert[tris[i][0]][2];
            if (vert[tris[i][0]][2] > max) max = vert[tris[i][0]][2];
        }

        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        // triangles
        glBegin(GL_TRIANGLES);
        for (int i = 0; i < m_scene->sceneSolution()->ordView().get_num_triangles(); i++)
        {
            int color = vert[tris[i][0]][2];
            glColor3d(paletteOrder[color][0], paletteOrder[color][1], paletteOrder[color][2]);

            glVertex2d(vert[tris[i][0]][0], vert[tris[i][0]][1]);
            glVertex2d(vert[tris[i][1]][0], vert[tris[i][1]][1]);
            glVertex2d(vert[tris[i][2]][0], vert[tris[i][2]][1]);
        }
        glEnd();

        glDisable(GL_POLYGON_OFFSET_FILL);

        glEndList();

        glCallList(m_listOrder);
    }
    else
    {
        glCallList(m_listOrder);
    }

    // paint labels
    if (Util::config()->orderLabel)
    {
        QFont fontLabel = font();
        fontLabel.setPointSize(fontLabel.pointSize() - 3);

        m_scene->sceneSolution()->ordView().lock_data();

        double3* vert = m_scene->sceneSolution()->ordView().get_vertices();
        int* lvert;
        char** ltext;
        double2* lbox;
        int nl = m_scene->sceneSolution()->ordView().get_labels(lvert, ltext, lbox);

        Point size((2.0/contextWidth()*fontMetrics().width(" "))/m_scale2d*aspect(),
                   (2.0/contextHeight()*fontMetrics().height())/m_scale2d);

        for (int i = 0; i < nl; i++)
        {
            glColor3f(1, 1, 1);
            // if (lbox[i][0]/m_scale*aspect() > size.x && lbox[i][1]/m_scale > size.y)
            {
                renderText(vert[lvert[i]][0] - size.x / 2.0,
                           vert[lvert[i]][1] - size.y / 2.0,
                           0.0,
                           ltext[i],
                           fontLabel);
            }
        }

        m_scene->sceneSolution()->ordView().unlock_data();
    }
}

void SceneView::paintOrderColorBar()
{
    if (!m_isSolutionPrepared) return;

    // order scalar view
    m_scene->sceneSolution()->ordView().lock_data();

    double3* vert = m_scene->sceneSolution()->ordView().get_vertices();
    int3* tris = m_scene->sceneSolution()->ordView().get_triangles();

    int min = 11;
    int max = 1;
    for (int i = 0; i < m_scene->sceneSolution()->ordView().get_num_triangles(); i++)
    {
        if (vert[tris[i][0]][2] < min) min = vert[tris[i][0]][2];
        if (vert[tris[i][0]][2] > max) max = vert[tris[i][0]][2];
    }

    m_scene->sceneSolution()->ordView().unlock_data();

    // order color map
    loadProjection2d();

    glScaled(2.0 / contextWidth(), 2.0 / contextHeight(), 1.0);
    glTranslated(- contextWidth() / 2.0, -contextHeight() / 2.0, 0.0);

    // dimensions
    int textWidth = fontMetrics().width("00");
    int textHeight = fontMetrics().height();
    Point scaleSize = Point(20 + 3 * textWidth, (20 + max * (2 * textHeight) - textHeight / 2.0 + 2));
    Point scaleBorder = Point(10.0, 10.0);
    double scaleLeft = (contextWidth() - (20 + 3 * textWidth));

    // blended rectangle
    drawBlend(Point(scaleLeft, scaleBorder.y), Point(scaleLeft + scaleSize.x - scaleBorder.x, scaleBorder.y + scaleSize.y),
              0.91, 0.91, 0.91);

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // bars
    glBegin(GL_QUADS);
    for (int i = 1; i < max+1; i++)
    {
        glColor3f(0.0, 0.0, 0.0);
        glVertex2d(scaleLeft + 10,                                 scaleBorder.y + 10 + (i-1)*(2 * textHeight));
        glVertex2d(scaleLeft + 10 + 3 * textWidth - scaleBorder.x, scaleBorder.y + 10 + (i-1)*(2 * textHeight));
        glVertex2d(scaleLeft + 10 + 3 * textWidth - scaleBorder.x, scaleBorder.y + 12 + (i )*(2 * textHeight) - textHeight / 2.0);
        glVertex2d(scaleLeft + 10,                                 scaleBorder.y + 12 + (i )*(2 * textHeight) - textHeight / 2.0);

        glColor3d(paletteOrder[i][0], paletteOrder[i][1], paletteOrder[i][2]);
        glVertex2d(scaleLeft + 12,                                     scaleBorder.y + 12 + (i-1)*(2 * textHeight));
        glVertex2d(scaleLeft + 10 + 3 * textWidth - 2 - scaleBorder.x, scaleBorder.y + 12 + (i-1)*(2 * textHeight));
        glVertex2d(scaleLeft + 10 + 3 * textWidth - 2 - scaleBorder.x, scaleBorder.y + 10 + (i  )*(2 * textHeight) - textHeight / 2.0);
        glVertex2d(scaleLeft + 12,                                     scaleBorder.y + 10 + (i  )*(2 * textHeight) - textHeight / 2.0);
    }
    glEnd();

    glDisable(GL_POLYGON_OFFSET_FILL);

    // labels
    glColor3f(1.0, 1.0, 1.0);
    for (int i = 1; i < max + 1; i++)
    {
        int sizeNumber = fontMetrics().width(QString::number(i));
        renderText(scaleLeft + 10 + 1.5 * textWidth - sizeNumber,
                   scaleBorder.y + 10.0 + (i-1)*(2.0 * textHeight) + textHeight / 2.0,
                   0.0,
                   QString::number(i));
    }
}

void SceneView::paintScalarFieldColorBar(double min, double max)
{
    loadProjection2d();

    glScaled(2.0 / contextWidth(), 2.0 / contextHeight(), 1.0);
    glTranslated(-contextWidth() / 2.0, -contextHeight() / 2.0, 0.0);

    // dimensions
    int textWidth = fontMetrics().width(QString::number(-1.0, '+e', 1)) + 3;
    int textHeight = fontMetrics().height();
    Point scaleSize = Point(45.0 + textWidth, contextHeight() - 20.0);
    Point scaleBorder = Point(10.0, 10.0);
    double scaleLeft = (contextWidth() - (45.0 + textWidth));
    int numTicks = 11;

    // blended rectangle
    drawBlend(Point(scaleLeft, scaleBorder.y), Point(scaleLeft + scaleSize.x - scaleBorder.x, scaleBorder.y + scaleSize.y),
              0.91, 0.91, 0.91);

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // palette border
    glColor3f(0.0, 0.0, 0.0);
    glBegin(GL_QUADS);
    glVertex2d(scaleLeft + 30.0, scaleBorder.y + scaleSize.y - 50.0);
    glVertex2d(scaleLeft + 10.0, scaleBorder.y + scaleSize.y - 50.0);
    glVertex2d(scaleLeft + 10.0, scaleBorder.y + 10.0);
    glVertex2d(scaleLeft + 30.0, scaleBorder.y + 10.0);
    glEnd();

    glDisable(GL_POLYGON_OFFSET_FILL);

    // palette
    glEnable(GL_TEXTURE_1D);
    glBindTexture(GL_TEXTURE_1D, 1);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

    glBegin(GL_QUADS);
    if (fabs(m_sceneViewSettings.scalarRangeMin - m_sceneViewSettings.scalarRangeMax) > EPS_ZERO)
        glTexCoord1d(m_texScale + m_texShift);
    else
        glTexCoord1d(m_texShift);
    glVertex2d(scaleLeft + 28.0, scaleBorder.y + scaleSize.y - 52.0);
    glVertex2d(scaleLeft + 12.0, scaleBorder.y + scaleSize.y - 52.0);
    glTexCoord1d(m_texShift);
    glVertex2d(scaleLeft + 12.0, scaleBorder.y + 12.0);
    glVertex2d(scaleLeft + 28.0, scaleBorder.y + 12.0);
    glEnd();

    glDisable(GL_TEXTURE_1D);

    // ticks
    glLineWidth(1.0);
    glBegin(GL_LINES);
    for (int i = 1; i < numTicks+1; i++)
    {
        double tickY = (scaleSize.y - 60.0) / (numTicks - 1.0);

        glVertex2d(scaleLeft + 10.0, scaleBorder.y + scaleSize.y - 49.0 - i*tickY);
        glVertex2d(scaleLeft + 15.0, scaleBorder.y + scaleSize.y - 49.0 - i*tickY);
        glVertex2d(scaleLeft + 25.0, scaleBorder.y + scaleSize.y - 49.0 - i*tickY);
        glVertex2d(scaleLeft + 30.0, scaleBorder.y + scaleSize.y - 49.0 - i*tickY);
    }
    glEnd();

    // labels
    for (int i = 1; i < numTicks+1; i++)
    {
        double value = 0.0;
        if (!Util::config()->scalarRangeLog)
            value = min + (double) (i-1) / (numTicks-1) * (max - min);
        else
            value = min + (double) pow(Util::config()->scalarRangeBase, ((i-1) / (numTicks-1)))/Util::config()->scalarRangeBase * (max - min);

        if (fabs(value) < EPS_ZERO) value = 0.0;
        double tickY = (scaleSize.y - 60.0) / (numTicks - 1.0);

        renderText(scaleLeft + 33.0 + ((value >= 0.0) ? fontMetrics().width("-") : 0.0),
                   scaleBorder.y + 10.0 + (i-1)*tickY - textHeight / 4.0,
                   0.0,
                   QString::number(value, '+e', 1));
    }

    // variable
    QString str = QString("%1 (%2)").
                  arg(physicFieldVariableShortcutString(m_sceneViewSettings.scalarPhysicFieldVariable)).
                  arg(physicFieldVariableUnitsString(m_sceneViewSettings.scalarPhysicFieldVariable));

    renderText(scaleLeft + scaleSize.x / 2.0 - fontMetrics().width(str) / 2.0,
               scaleBorder.y + scaleSize.y - 20.0,
               0.0,
               str);
    // line
    glLineWidth(1.0);
    glBegin(GL_LINES);
    glVertex2d(scaleLeft + 5.0, scaleBorder.y + scaleSize.y - 31.0);
    glVertex2d(scaleLeft + scaleSize.x - 15.0, scaleBorder.y + scaleSize.y - 31.0);
    glEnd();
}

void SceneView::paintScalarField()
{
    if (!m_isSolutionPrepared) return;

    loadProjection2d(true);

    if (m_listScalarField == -1)
    {
        qDebug() << "SceneView::paintScalarField(), min = " << m_sceneViewSettings.scalarRangeMin << ", max = " << m_sceneViewSettings.scalarRangeMax;

        m_listScalarField = glGenLists(1);
        glNewList(m_listScalarField, GL_COMPILE);

        // range
        double irange = 1.0 / (m_sceneViewSettings.scalarRangeMax - m_sceneViewSettings.scalarRangeMin);
        // special case: constant solution
        if (fabs(m_sceneViewSettings.scalarRangeMax - m_sceneViewSettings.scalarRangeMin) < EPS_ZERO)
            irange = 1.0;

        m_scene->sceneSolution()->linScalarView().lock_data();

        double3* linVert = m_scene->sceneSolution()->linScalarView().get_vertices();
        int3* linTris = m_scene->sceneSolution()->linScalarView().get_triangles();
        Point point[3];
        double value[3];

        // set texture for coloring
        glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D, 1);
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

        // set texture transformation matrix
        glMatrixMode(GL_TEXTURE);
        glLoadIdentity();
        glTranslated(m_texShift, 0.0, 0.0);
        glScaled(m_texScale, 0.0, 0.0);

        glBegin(GL_TRIANGLES);
        for (int i = 0; i < m_scene->sceneSolution()->linScalarView().get_num_triangles(); i++)
        {
            for (int j = 0; j < 3; j++)
            {
                point[j].x = linVert[linTris[i][j]][0];
                point[j].y = linVert[linTris[i][j]][1];
                value[j]   = linVert[linTris[i][j]][2];
            }

            if (!m_sceneViewSettings.scalarRangeAuto)
            {
                double avgValue = (value[0] + value[1] + value[2]) / 3.0;
                if (avgValue < m_sceneViewSettings.scalarRangeMin || avgValue > m_sceneViewSettings.scalarRangeMax)
                    continue;
            }

            for (int j = 0; j < 3; j++)
            {
                if (Util::config()->scalarRangeLog)
                    glTexCoord1d(log10(1.0 + (Util::config()->scalarRangeBase-1.0)*(value[j] - m_sceneViewSettings.scalarRangeMin) * irange)/log10(Util::config()->scalarRangeBase));
                else
                    glTexCoord1d((value[j] - m_sceneViewSettings.scalarRangeMin) * irange);
                glVertex2d(point[j].x, point[j].y);
            }
        }
        glEnd();

        glDisable(GL_POLYGON_OFFSET_FILL);
        glDisable(GL_TEXTURE_1D);

        // switch-off texture transform
        glMatrixMode(GL_TEXTURE);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);

        m_scene->sceneSolution()->linScalarView().unlock_data();

        glEndList();

        glCallList(m_listScalarField);
    }
    else
    {
        glCallList(m_listScalarField);
    }
}

void SceneView::paintScalarField3D()
{
    if (!m_isSolutionPrepared) return;

    loadProjection3d(true);

    if (m_listScalarField3D == -1)
    {
        qDebug() << "SceneView::paintScalarField3D(), min = " << m_sceneViewSettings.scalarRangeMin << ", max = " << m_sceneViewSettings.scalarRangeMax;

        m_listScalarField3D = glGenLists(1);
        glNewList(m_listScalarField3D, GL_COMPILE);

        // gradient background
        paintBackground();
        glEnable(GL_DEPTH_TEST);

        // range
        double irange = 1.0 / (m_sceneViewSettings.scalarRangeMax - m_sceneViewSettings.scalarRangeMin);
        // special case: constant solution
        if (fabs(m_sceneViewSettings.scalarRangeMin - m_sceneViewSettings.scalarRangeMax) < EPS_ZERO)
        {
            irange = 1.0;
        }

        m_scene->sceneSolution()->linScalarView().lock_data();

        double3* linVert = m_scene->sceneSolution()->linScalarView().get_vertices();
        int3* linTris = m_scene->sceneSolution()->linScalarView().get_triangles();
        Point point[3];
        double value[3];

        double max = qMax(m_scene->boundingBox().width(), m_scene->boundingBox().height());

        if (!Util::config()->scalarView3DLighting)
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

        glPushMatrix();
        glScaled(1.0, 1.0, max/4.0 * 1.0/(fabs(m_sceneViewSettings.scalarRangeMin - m_sceneViewSettings.scalarRangeMax)));

        initLighting();

        // set texture for coloring
        glEnable(GL_TEXTURE_1D);
        glBindTexture(GL_TEXTURE_1D, 1);

        // set texture transformation matrix
        glMatrixMode(GL_TEXTURE);
        glLoadIdentity();
        glTranslated(m_texShift, 0.0, 0.0);
        glScaled(m_texScale, 0.0, 0.0);

        glBegin(GL_TRIANGLES);
        for (int i = 0; i < m_scene->sceneSolution()->linScalarView().get_num_triangles(); i++)
        {
            point[0].x = linVert[linTris[i][0]][0];
            point[0].y = linVert[linTris[i][0]][1];
            value[0]   = linVert[linTris[i][0]][2];
            point[1].x = linVert[linTris[i][1]][0];
            point[1].y = linVert[linTris[i][1]][1];
            value[1]   = linVert[linTris[i][1]][2];
            point[2].x = linVert[linTris[i][2]][0];
            point[2].y = linVert[linTris[i][2]][1];
            value[2]   = linVert[linTris[i][2]][2];

            if (!m_sceneViewSettings.scalarRangeAuto)
            {
                double avgValue = (value[0] + value[1] + value[2]) / 3.0;
                if (avgValue < m_sceneViewSettings.scalarRangeMin || avgValue > m_sceneViewSettings.scalarRangeMax)
                    continue;
            }

            double delta = 0.0;

            if (Util::config()->scalarView3DLighting)
            {
                double* normal = computeNormal(point[0].x, point[0].y, - delta - (value[0] - m_sceneViewSettings.scalarRangeMin),
                                               point[1].x, point[1].y, - delta - (value[1] - m_sceneViewSettings.scalarRangeMin),
                                               point[2].x, point[2].y, - delta - (value[2] - m_sceneViewSettings.scalarRangeMin));
                glNormal3d(normal[0], normal[1], normal[2]);
            }
            for (int j = 0; j < 3; j++)
            {
                glTexCoord1d((value[j] - m_sceneViewSettings.scalarRangeMin) * irange);
                glVertex3d(point[j].x, point[j].y, - delta - (value[j] - m_sceneViewSettings.scalarRangeMin));
            }
        }
        glEnd();

        glDisable(GL_TEXTURE_1D);
        glDisable(GL_LIGHTING);

        // draw blended mesh
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glColor4d(0.5, 0.5, 0.5, 0.3);

        glBegin(GL_TRIANGLES);
        for (int i = 0; i < m_scene->sceneSolution()->meshInitial()->get_num_active_elements(); i++)
        {
            Element *element = m_scene->sceneSolution()->meshInitial()->get_element(i);
            if (element->is_triangle())
            {
                glVertex3d(element->vn[0]->x, element->vn[0]->y, 0.0);
                glVertex3d(element->vn[1]->x, element->vn[1]->y, 0.0);
                glVertex3d(element->vn[2]->x, element->vn[2]->y, 0.0);
            }
        }
        glEnd();

        glDisable(GL_POLYGON_OFFSET_FILL);
        glDisable(GL_BLEND);

        // geometry - edges
        foreach (SceneEdge *edge, m_scene->edges)
        {

            glColor3f(Util::config()->colorEdges.redF(),
                      Util::config()->colorEdges.greenF(),
                      Util::config()->colorEdges.blueF());
            glLineWidth(Util::config()->edgeWidth);

            if (edge->angle == 0)
            {
                glBegin(GL_LINES);
                glVertex3d(edge->nodeStart->point.x, edge->nodeStart->point.y, 0.0);
                glVertex3d(edge->nodeEnd->point.x, edge->nodeEnd->point.y, 0.0);
                glEnd();
            }
            else
            {
                Point center = edge->center();
                double radius = edge->radius();
                double startAngle = atan2(center.y - edge->nodeStart->point.y, center.x - edge->nodeStart->point.x) / M_PI*180.0 - 180.0;

                drawArc(center, radius, startAngle, edge->angle, edge->angle/2);
            }

            glDisable(GL_LINE_STIPPLE);
            glLineWidth(1.0);
        }

        glDisable(GL_DEPTH_TEST);

        // switch-off texture transform
        glMatrixMode(GL_TEXTURE);
        glLoadIdentity();
        glMatrixMode(GL_MODELVIEW);

        glPopMatrix();

        m_scene->sceneSolution()->linScalarView().unlock_data();

        glEndList();

        glCallList(m_listScalarField3D);
    }
    else
    {
        glCallList(m_listScalarField3D);
    }
}

void SceneView::paintScalarField3DSolid()
{
    if (!m_isSolutionPrepared) return;

    loadProjection3d(true);

    if (m_listScalarField3DSolid == -1)
    {
        qDebug() << "SceneView::paintScalarField3DSolid(), min = " << m_sceneViewSettings.scalarRangeMin << ", max = " << m_sceneViewSettings.scalarRangeMax;

        m_listScalarField3DSolid = glGenLists(1);
        glNewList(m_listScalarField3DSolid, GL_COMPILE);

        bool isModel = (m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_Model);

        // gradient background
        paintBackground();
        glEnable(GL_DEPTH_TEST);

        RectPoint rect = m_scene->boundingBox();
        double max = qMax(rect.width(), rect.height());
        double depth = max / 4.0;

        // range
        double irange = 1.0 / (m_sceneViewSettings.scalarRangeMax - m_sceneViewSettings.scalarRangeMin);
        // special case: constant solution
        if (fabs(m_sceneViewSettings.scalarRangeMin - m_sceneViewSettings.scalarRangeMax) < EPS_ZERO)
        {
            irange = 1.0;
        }

        double phi = Util::config()->scalarView3DAngle;

        m_scene->sceneSolution()->linScalarView().lock_data();

        double3* linVert = m_scene->sceneSolution()->linScalarView().get_vertices();
        int3* linTris = m_scene->sceneSolution()->linScalarView().get_triangles();
        int3* linEdges = m_scene->sceneSolution()->linScalarView().get_edges();
        Point point[3];
        double value[3];
        double *normal;

        glPushMatrix();

        // set texture for coloring
        if (!isModel)
        {
            glEnable(GL_TEXTURE_1D);
            glBindTexture(GL_TEXTURE_1D, 1);
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);

            // set texture transformation matrix
            glMatrixMode(GL_TEXTURE);
            glLoadIdentity();
            glTranslated(m_texShift, 0.0, 0.0);
            glScaled(m_texScale, 0.0, 0.0);
        }

        initLighting();

        if (m_scene->problemInfo()->problemType == ProblemType_Planar)
        {
            glBegin(GL_TRIANGLES);
            for (int i = 0; i < m_scene->sceneSolution()->linScalarView().get_num_triangles(); i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    point[j].x = linVert[linTris[i][j]][0];
                    point[j].y = linVert[linTris[i][j]][1];
                    value[j]   = linVert[linTris[i][j]][2];
                }

                if (!m_sceneViewSettings.scalarRangeAuto)
                {
                    double avgValue = (value[0] + value[1] + value[2]) / 3.0;
                    if (avgValue < m_sceneViewSettings.scalarRangeMin || avgValue > m_sceneViewSettings.scalarRangeMax)
                        continue;
                }

                // z = - depth / 2.0
                normal = computeNormal(point[0].x, point[0].y, -depth/2.0, point[1].x, point[1].y, -depth/2.0, point[2].x, point[2].y, -depth/2.0);
                glNormal3d(normal[0], normal[1], normal[2]);

                for (int j = 0; j < 3; j++)
                {
                    if (!isModel) glTexCoord1d((value[j] - m_sceneViewSettings.scalarRangeMin) * irange);
                    glVertex3d(point[j].x, point[j].y, -depth/2.0);
                }

                // z = + depth / 2.0
                normal = computeNormal(point[0].x, point[0].y, depth/2.0, point[1].x, point[1].y, depth/2.0, point[2].x, point[2].y, depth/2.0);
                glNormal3d(normal[0], normal[1], normal[2]);

                for (int j = 0; j < 3; j++)
                {
                    if (!isModel) glTexCoord1d((value[j] - m_sceneViewSettings.scalarRangeMin) * irange);
                    glVertex3d(point[j].x, point[j].y, depth/2.0);
                }
            }
            glEnd();

            // length
            glBegin(GL_QUADS);
            for (int i = 0; i < m_scene->sceneSolution()->linScalarView().get_num_edges(); i++)
            {
                // draw only boundary edges
                if (!linEdges[i][2]) continue;

                for (int j = 0; j < 2; j++)
                {
                    point[j].x = linVert[linEdges[i][j]][0];
                    point[j].y = linVert[linEdges[i][j]][1];
                    value[j]   = linVert[linEdges[i][j]][2];
                }

                if (!m_sceneViewSettings.scalarRangeAuto)
                {
                    double avgValue = (value[0] + value[1] + value[2]) / 3.0;
                    if (avgValue < m_sceneViewSettings.scalarRangeMin || avgValue > m_sceneViewSettings.scalarRangeMax)
                        continue;
                }

                normal = computeNormal(point[0].x, point[0].y, -depth/2.0, point[1].x, point[1].y, -depth/2.0, point[1].x, point[1].y,  depth/2.0);
                glNormal3d(normal[0], normal[1], normal[2]);

                if (!isModel) glTexCoord1d((value[0] - m_sceneViewSettings.scalarRangeMin) * irange);
                glVertex3d(point[0].x, point[0].y, -depth/2.0);
                if (!isModel) glTexCoord1d((value[1] - m_sceneViewSettings.scalarRangeMin) * irange);
                glVertex3d(point[1].x, point[1].y, -depth/2.0);
                if (!isModel) glTexCoord1d((value[1] - m_sceneViewSettings.scalarRangeMin) * irange);
                glVertex3d(point[1].x, point[1].y, depth/2.0);
                if (!isModel) glTexCoord1d((value[0] - m_sceneViewSettings.scalarRangeMin) * irange);
                glVertex3d(point[0].x, point[0].y, depth/2.0);
            }
            glEnd();
        }
        else
        {
            // side
            glBegin(GL_TRIANGLES);
            for (int i = 0; i < m_scene->sceneSolution()->linScalarView().get_num_triangles(); i++)
            {
                for (int j = 0; j < 3; j++)
                {
                    point[j].x = linVert[linTris[i][j]][0];
                    point[j].y = linVert[linTris[i][j]][1];
                    value[j]   = linVert[linTris[i][j]][2];
                }

                if (!m_sceneViewSettings.scalarRangeAuto)
                {
                    double avgValue = (value[0] + value[1] + value[2]) / 3.0;
                    if (avgValue < m_sceneViewSettings.scalarRangeMin || avgValue > m_sceneViewSettings.scalarRangeMax)
                        continue;
                }

                for (int j = 0; j < 2; j++)
                {
                    normal = computeNormal(point[0].x * cos(j*phi/180.0*M_PI), point[0].y, point[0].x * sin(j*phi/180.0*M_PI),
                                           point[1].x * cos(j*phi/180.0*M_PI), point[1].y, point[1].x * sin(j*phi/180.0*M_PI),
                                           point[2].x * cos(j*phi/180.0*M_PI), point[2].y, point[2].x * sin(j*phi/180.0*M_PI));
                    glNormal3d(normal[0], normal[1], normal[2]);

                    glTexCoord1d((value[0] - m_sceneViewSettings.scalarRangeMin) * irange);
                    glVertex3d(point[0].x * cos(j*phi/180.0*M_PI), point[0].y, point[0].x * sin(j*phi/180.0*M_PI));
                    glTexCoord1d((value[1] - m_sceneViewSettings.scalarRangeMin) * irange);
                    glVertex3d(point[1].x * cos(j*phi/180.0*M_PI), point[1].y, point[1].x * sin(j*phi/180.0*M_PI));
                    glTexCoord1d((value[2] - m_sceneViewSettings.scalarRangeMin) * irange);
                    glVertex3d(point[2].x * cos(j*phi/180.0*M_PI), point[2].y, point[2].x * sin(j*phi/180.0*M_PI));
                }
            }
            glEnd();

            // symmetry
            glBegin(GL_QUADS);
            for (int i = 0; i < m_scene->sceneSolution()->linScalarView().get_num_edges(); i++)
            {
                // draw only boundary edges
                if (!linEdges[i][2]) continue;

                for (int j = 0; j < 2; j++)
                {
                    point[j].x = linVert[linEdges[i][j]][0];
                    point[j].y = linVert[linEdges[i][j]][1];
                    value[j]   = linVert[linEdges[i][j]][2];
                }

                if (!m_sceneViewSettings.scalarRangeAuto)
                {
                    double avgValue = (value[0] + value[1] + value[2]) / 3.0;
                    if (avgValue < m_sceneViewSettings.scalarRangeMin || avgValue > m_sceneViewSettings.scalarRangeMax)
                        continue;
                }

                int count = 30;
                double step = phi/count;
                for (int j = 0; j < count; j++)
                {
                    normal = computeNormal(point[0].x * cos((j+0)*step/180.0*M_PI), point[0].y, point[0].x * sin((j+0)*step/180.0*M_PI),
                                           point[1].x * cos((j+0)*step/180.0*M_PI), point[1].y, point[1].x * sin((j+0)*step/180.0*M_PI),
                                           point[1].x * cos((j+1)*step/180.0*M_PI), point[1].y, point[1].x * sin((j+1)*step/180.0*M_PI));
                    glNormal3d(normal[0], normal[1], normal[2]);

                    if (!isModel) glTexCoord1d((value[0] - m_sceneViewSettings.scalarRangeMin) * irange);
                    glVertex3d(point[0].x * cos((j+0)*step/180.0*M_PI), point[0].y, point[0].x * sin((j+0)*step/180.0*M_PI));
                    if (!isModel) glTexCoord1d((value[1] - m_sceneViewSettings.scalarRangeMin) * irange);
                    glVertex3d(point[1].x * cos((j+0)*step/180.0*M_PI), point[1].y, point[1].x * sin((j+0)*step/180.0*M_PI));
                    if (!isModel) glTexCoord1d((value[1] - m_sceneViewSettings.scalarRangeMin) * irange);
                    glVertex3d(point[1].x * cos((j+1)*step/180.0*M_PI), point[1].y, point[1].x * sin((j+1)*step/180.0*M_PI));
                    if (!isModel) glTexCoord1d((value[0] - m_sceneViewSettings.scalarRangeMin) * irange);
                    glVertex3d(point[0].x * cos((j+1)*step/180.0*M_PI), point[0].y, point[0].x * sin((j+1)*step/180.0*M_PI));
                }
            }
            glEnd();
        }

        glDisable(GL_POLYGON_OFFSET_FILL);
        glDisable(GL_LIGHTING);

        if (!isModel)
        {
            glDisable(GL_TEXTURE_1D);

            // switch-off texture transform
            glMatrixMode(GL_TEXTURE);
            glLoadIdentity();
            glMatrixMode(GL_MODELVIEW);
        }

        // geometry
        if (m_scene->problemInfo()->problemType == ProblemType_Planar)
        {
            glColor3f(Util::config()->colorEdges.redF(),
                      Util::config()->colorEdges.greenF(),
                      Util::config()->colorEdges.blueF());
            glLineWidth(Util::config()->edgeWidth);

            foreach (SceneEdge *edge, m_scene->edges)
            {
                for (int j = 0; j < 2; j++)
                {
                    if (edge->angle == 0)
                    {
                        glBegin(GL_LINES);
                        glVertex3d(edge->nodeStart->point.x, edge->nodeStart->point.y, - depth/2.0 + j*depth);
                        glVertex3d(edge->nodeEnd->point.x, edge->nodeEnd->point.y, - depth/2.0 + j*depth);
                        glEnd();
                    }
                    else
                    {
                        Point center = edge->center();
                        double radius = edge->radius();
                        double startAngle = atan2(center.y - edge->nodeStart->point.y, center.x - edge->nodeStart->point.x) / M_PI*180.0 - 180.0;

                        double theta = edge->angle / double(edge->angle/2 - 1);

                        glBegin(GL_LINE_STRIP);
                        for (int i = 0; i < edge->angle/2; i++)
                        {
                            double arc = (startAngle + i*theta)/180.0*M_PI;

                            double x = radius * cos(arc);
                            double y = radius * sin(arc);

                            glVertex3d(center.x + x, center.y + y, - depth/2.0 + j*depth);
                        }
                        glEnd();
                    }
                }
            }
            glLineWidth(1.0);
        }
        else
        {
            // geometry
            glColor3d(Util::config()->colorEdges.redF(),
                      Util::config()->colorEdges.greenF(),
                      Util::config()->colorEdges.blueF());
            glLineWidth(Util::config()->edgeWidth);

            foreach (SceneEdge *edge, m_scene->edges)
            {
                for (int j = 0; j < 2; j++)
                {
                    if (edge->angle == 0)
                    {
                        glBegin(GL_LINES);
                        glVertex3d(edge->nodeStart->point.x * cos(j*phi/180.0*M_PI), edge->nodeStart->point.y, edge->nodeStart->point.x * sin(j*phi/180.0*M_PI));
                        glVertex3d(edge->nodeEnd->point.x * cos(j*phi/180.0*M_PI), edge->nodeEnd->point.y, edge->nodeEnd->point.x * sin(j*phi/180.0*M_PI));
                        glEnd();
                    }
                    else
                    {
                        Point center = edge->center();
                        double radius = edge->radius();
                        double startAngle = atan2(center.y - edge->nodeStart->point.y, center.x - edge->nodeStart->point.x) / M_PI*180.0 - 180.0;

                        double theta = edge->angle / double(edge->angle/2 - 1);

                        glBegin(GL_LINE_STRIP);
                        for (int i = 0; i < edge->angle/2; i++)
                        {
                            double arc = (startAngle + i*theta)/180.0*M_PI;

                            double x = radius * cos(arc);
                            double y = radius * sin(arc);

                            glVertex3d((center.x + x) * cos(j*phi/180.0*M_PI), center.y + y, (center.x + x) * sin(j*phi/180.0*M_PI));
                        }
                        glEnd();
                    }
                }
            }
            glLineWidth(1.0);
        }

        glDisable(GL_DEPTH_TEST);

        glPopMatrix();

        m_scene->sceneSolution()->linScalarView().unlock_data();

        glEndList();

        glCallList(m_listScalarField3DSolid);
    }
    else
    {
        glCallList(m_listScalarField3DSolid);
    }
}

void SceneView::paintContours()
{
    if (!m_isSolutionPrepared) return;

    loadProjection2d(true);

    if (m_listContours == -1)
    {
        m_listContours = glGenLists(1);
        glNewList(m_listContours, GL_COMPILE);

        m_scene->sceneSolution()->linContourView().lock_data();

        double3* tvert = m_scene->sceneSolution()->linContourView().get_vertices();
        int3* tris = m_scene->sceneSolution()->linContourView().get_triangles();

        // transform variable
        double rangeMin =  CONST_DOUBLE;
        double rangeMax = -CONST_DOUBLE;

        double3* vert = new double3[m_scene->sceneSolution()->linContourView().get_num_vertices()];
        for (int i = 0; i < m_scene->sceneSolution()->linContourView().get_num_vertices(); i++)
        {
            vert[i][0] = tvert[i][0];
            vert[i][1] = tvert[i][1];
            vert[i][2] = tvert[i][2];

            if (vert[i][2] > rangeMax) rangeMax = tvert[i][2];
            if (vert[i][2] < rangeMin) rangeMin = tvert[i][2];
        }

        qDebug() << "SceneView::paintContours(), min = " << rangeMin << ", max = " << rangeMax;

        // value range
        double step = (rangeMax-rangeMin)/Util::config()->contoursCount;

        // draw contours
        glLineWidth(1.0);
        glColor3f(Util::config()->colorContours.redF(),
                  Util::config()->colorContours.greenF(),
                  Util::config()->colorContours.blueF());

        glBegin(GL_LINES);
        for (int i = 0; i < m_scene->sceneSolution()->linContourView().get_num_triangles(); i++)
        {
            if (finite(vert[tris[i][0]][2]) && finite(vert[tris[i][1]][2]) && finite(vert[tris[i][2]][2]))
            {
                paintContoursTri(vert, &tris[i], step);
            }
        }
        glEnd();

        delete vert;

        m_scene->sceneSolution()->linContourView().unlock_data();

        glEndList();

        glCallList(m_listContours);
    }
    else
    {
        glCallList(m_listContours);
    }
}

void SceneView::paintContoursTri(double3* vert, int3* tri, double step)
{
    // sort the vertices by their value, keep track of the permutation sign
    int i, idx[3], perm = 0;
    memcpy(idx, tri, sizeof(idx));
    for (i = 0; i < 2; i++)
    {
        if (vert[idx[0]][2] > vert[idx[1]][2]) { std::swap(idx[0], idx[1]); perm++; }
        if (vert[idx[1]][2] > vert[idx[2]][2]) { std::swap(idx[1], idx[2]); perm++; }
    }
    if (fabs(vert[idx[0]][2] - vert[idx[2]][2]) < 1e-3 * fabs(step)) return;

    // get the first (lowest) contour value
    double val = vert[idx[0]][2];

    double y = ceil(val / step);
    if (y < val / step) y + 1.0;
    val = y * step;

    int l1 = 0, l2 = 1;
    int r1 = 0, r2 = 2;
    while (val < vert[idx[r2]][2])
    {
        double ld = vert[idx[l2]][2] - vert[idx[l1]][2];
        double rd = vert[idx[r2]][2] - vert[idx[r1]][2];

        // draw a slice of the triangle
        while (val < vert[idx[l2]][2])
        {
            double lt = (val - vert[idx[l1]][2]) / ld;
            double rt = (val - vert[idx[r1]][2]) / rd;

            double x1 = (1.0 - lt) * vert[idx[l1]][0] + lt * vert[idx[l2]][0];
            double y1 = (1.0 - lt) * vert[idx[l1]][1] + lt * vert[idx[l2]][1];
            double x2 = (1.0 - rt) * vert[idx[r1]][0] + rt * vert[idx[r2]][0];
            double y2 = (1.0 - rt) * vert[idx[r1]][1] + rt * vert[idx[r2]][1];

            if (perm & 1) { glVertex2d(x1, y1); glVertex2d(x2, y2); }
            else { glVertex2d(x2, y2); glVertex2d(x1, y1); }

            val += step;
        }
        l1 = 1;
        l2 = 2;
    }
}

void SceneView::paintVectors()
{
    if (!m_isSolutionPrepared) return;

    loadProjection2d(true);

    if (m_listVectors == -1)
    {
        m_listVectors = glGenLists(1);
        glNewList(m_listVectors, GL_COMPILE);

        double vectorRangeMin = m_scene->sceneSolution()->vecVectorView().get_min_value();
        double vectorRangeMax = m_scene->sceneSolution()->vecVectorView().get_max_value();

        qDebug() << "SceneView::paintVectors(), min = " << vectorRangeMin << ", max = " << vectorRangeMax;

        double irange = 1.0 / (vectorRangeMax - vectorRangeMin);
        if (fabs(vectorRangeMin - vectorRangeMax) < EPS_ZERO) return;

        RectPoint rect = m_scene->boundingBox();
        double gs = (rect.width() + rect.height()) / Util::config()->vectorCount;

        // paint
        m_scene->sceneSolution()->vecVectorView().lock_data();

        double4* vecVert = m_scene->sceneSolution()->vecVectorView().get_vertices();
        int3* vecTris = m_scene->sceneSolution()->vecVectorView().get_triangles();

        glEnable(GL_POLYGON_OFFSET_FILL);
        glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

        glBegin(GL_TRIANGLES);
        for (int i = 0; i < m_scene->sceneSolution()->vecVectorView().get_num_triangles(); i++)
        {
            Point a(vecVert[vecTris[i][0]][0], vecVert[vecTris[i][0]][1]);
            Point b(vecVert[vecTris[i][1]][0], vecVert[vecTris[i][1]][1]);
            Point c(vecVert[vecTris[i][2]][0], vecVert[vecTris[i][2]][1]);

            RectPoint r;
            r.start = Point(qMin(qMin(a.x, b.x), c.x), qMin(qMin(a.y, b.y), c.y));
            r.end = Point(qMax(qMax(a.x, b.x), c.x), qMax(qMax(a.y, b.y), c.y));

            // double area
            double area2 = a.x * (b.y - c.y) + b.x * (c.y - a.y) + c.x * (a.y - b.y);

            // plane equation
            double aa = b.x*c.y - c.x*b.y;
            double ab = c.x*a.y - a.x*c.y;
            double ac = a.x*b.y - b.x*a.y;
            double ba = b.y - c.y;
            double bb = c.y - a.y;
            double bc = a.y - b.y;
            double ca = c.x - b.x;
            double cb = a.x - c.x;
            double cc = b.x - a.x;

            double ax = (aa * vecVert[vecTris[i][0]][2] + ab * vecVert[vecTris[i][1]][2] + ac * vecVert[vecTris[i][2]][2]) / area2;
            double bx = (ba * vecVert[vecTris[i][0]][2] + bb * vecVert[vecTris[i][1]][2] + bc * vecVert[vecTris[i][2]][2]) / area2;
            double cx = (ca * vecVert[vecTris[i][0]][2] + cb * vecVert[vecTris[i][1]][2] + cc * vecVert[vecTris[i][2]][2]) / area2;

            double ay = (aa * vecVert[vecTris[i][0]][3] + ab * vecVert[vecTris[i][1]][3] + ac * vecVert[vecTris[i][2]][3]) / area2;
            double by = (ba * vecVert[vecTris[i][0]][3] + bb * vecVert[vecTris[i][1]][3] + bc * vecVert[vecTris[i][2]][3]) / area2;
            double cy = (ca * vecVert[vecTris[i][0]][3] + cb * vecVert[vecTris[i][1]][3] + cc * vecVert[vecTris[i][2]][3]) / area2;

            for (int j = floor(r.start.x / gs); j < ceil(r.end.x / gs); j++)
            {
                for (int k = floor(r.start.y / gs); k < ceil(r.end.y / gs); k++)
                {
                    Point point(j*gs, k*gs);
                    if (k % 2 == 0) point.x += gs/2.0;

                    // find in triangle
                    bool inTriangle = true;

                    for (int l = 0; l < 3; l++)
                    {
                        int p = l + 1;
                        if (p == 3)
                            p = 0;

                        double z = (vecVert[vecTris[i][p]][0] - vecVert[vecTris[i][l]][0]) * (point.y - vecVert[vecTris[i][l]][1]) -
                                   (vecVert[vecTris[i][p]][1] - vecVert[vecTris[i][l]][1]) * (point.x - vecVert[vecTris[i][l]][0]);

                        if (z < 0)
                        {
                            inTriangle = false;
                            break;
                        }
                    }

                    if (inTriangle)
                    {
                        // view
                        double dx = ax + bx * point.x + cx * point.y;
                        double dy = ay + by * point.x + cy * point.y;

                        double value = sqrt(sqr(dx) + sqr(dy));
                        double angle = atan2(dy, dx);

                        if (Util::config()->vectorProportional)
                        {
                            dx = ((value - vectorRangeMin) * irange) * Util::config()->vectorScale * gs * cos(angle);
                            dy = ((value - vectorRangeMin) * irange) * Util::config()->vectorScale * gs * sin(angle);
                        }
                        else
                        {
                            dx = Util::config()->vectorScale * gs * cos(angle);
                            dy = Util::config()->vectorScale * gs * sin(angle);
                        }

                        double dm = sqrt(sqr(dx) + sqr(dy));

                        // color
                        if (Util::config()->vectorColor)
                        {
                            double color = 0.7 - 0.7 * (value - vectorRangeMin) * irange;
                            glColor3f(color, color, color);
                        }
                        else
                        {
                            glColor3f(Util::config()->colorVectors.redF(),
                                      Util::config()->colorVectors.greenF(),
                                      Util::config()->colorVectors.blueF());
                        }

                        glVertex2d(point.x + dm/5.0 * cos(angle - M_PI_2), point.y + dm/5.0 * sin(angle - M_PI_2));
                        glVertex2d(point.x + dm/5.0 * cos(angle + M_PI_2), point.y + dm/5.0 * sin(angle + M_PI_2));
                        glVertex2d(point.x + dm     * cos(angle),          point.y + dm     * sin(angle));
                    }
                }
            }
        }
        glEnd();

        glDisable(GL_POLYGON_OFFSET_FILL);

        m_scene->sceneSolution()->vecVectorView().unlock_data();

        glEndList();

        glCallList(m_listVectors);
    }
    else
    {
        glCallList(m_listVectors);
    }
}

void SceneView::paintSceneModeLabel()
{
    QString text = "";

    switch (m_sceneMode)
    {
    case SceneMode_OperateOnNodes:
        text = tr("Operate on nodes");
        break;
    case SceneMode_OperateOnEdges:
        text = tr("Operate on edges");
        break;
    case SceneMode_OperateOnLabels:
        text = tr("Operate on labels");
        break;
    case SceneMode_Postprocessor:
        switch (m_sceneViewSettings.postprocessorShow)
        {
        case SceneViewPostprocessorShow_ScalarView:
        case SceneViewPostprocessorShow_ScalarView3D:
        case SceneViewPostprocessorShow_ScalarView3DSolid:
            text = physicFieldVariableString(m_sceneViewSettings.scalarPhysicFieldVariable);
            if (m_sceneViewSettings.scalarPhysicFieldVariableComp != PhysicFieldVariableComp_Scalar)
                text += " - " + physicFieldVariableCompString(m_sceneViewSettings.scalarPhysicFieldVariableComp);
            break;
        case SceneViewPostprocessorShow_Model:
            text = tr("Model");
            break;
        case SceneViewPostprocessorShow_Order:
            text = tr("Polynomial order");
            break;
        default:
            text = tr("Postprocessor");
        }
        break;
    }

    loadProjection2d();

    glLoadIdentity();

    glScaled(2.0/contextWidth(), 2.0/contextHeight(), 1.0);
    glTranslated(-contextWidth() / 2.0, -contextHeight() / 2.0, 0.0);

    glDisable(GL_DEPTH_TEST);

    // render viewport label
    QFont fontLabel = font();
    fontLabel.setPointSize(fontLabel.pointSize() + 1);

    Point posText = Point((contextWidth()-QFontMetrics(fontLabel).width(text)) / 2.0,
                          (contextHeight() - QFontMetrics(fontLabel).height() / 1.3));

    // blended rectangle
    double xs = posText.x - QFontMetrics(fontLabel).width(" ");
    double ys = posText.y - QFontMetrics(fontLabel).height() / 3.0;
    double xe = xs + QFontMetrics(fontLabel).width(text + "  ");
    double ye = contextHeight();

    drawBlend(Point(xs, ys), Point(xe, ye), 0.8, 0.8, 0.8, 0.93);

    // text
    glColor3f(0.0, 0.0, 0.0);
    renderText(posText.x, posText.y, 0.0, text, fontLabel);
}

void SceneView::paintZoomRegion()
{
    loadProjection2d(true);

    // zoom or select region
    if (m_region)
    {
        Point posStart = position(Point(m_regionPos.x(), m_regionPos.y()));
        Point posEnd = position(Point(m_lastPos.x(), m_lastPos.y()));

        drawBlend(posStart, posEnd,
                  Util::config()->colorHighlighted.redF(),
                  Util::config()->colorHighlighted.greenF(),
                  Util::config()->colorHighlighted.blueF());
    }
}

void SceneView::paintSnapToGrid()
{
    if (m_snapToGrid)
    {
        loadProjection2d(true);

        Point p = position(Point(m_lastPos.x(), m_lastPos.y()));

        Point snapPoint;
        snapPoint.x = round(p.x / Util::config()->gridStep) * Util::config()->gridStep;
        snapPoint.y = round(p.y / Util::config()->gridStep) * Util::config()->gridStep;

        glColor3f(Util::config()->colorHighlighted.redF(),
                  Util::config()->colorHighlighted.greenF(),
                  Util::config()->colorHighlighted.blueF());
        glPointSize(Util::config()->nodeSize - 1.0);
        glBegin(GL_POINTS);
        glVertex2d(snapPoint.x, snapPoint.y);
        glEnd();
    }
}

void SceneView::paintChartLine()
{
    loadProjection2d(true);

    glColor3f(Util::config()->colorSelected.redF(),
              Util::config()->colorSelected.greenF(),
              Util::config()->colorSelected.blueF());
    glLineWidth(3.0);

    glBegin(GL_LINES);
    glVertex2d(m_chartLine.start.x, m_chartLine.start.y);
    glVertex2d(m_chartLine.end.x, m_chartLine.end.y);
    glEnd();
}

const float* SceneView::paletteColor(double x)
{
    switch (Util::config()->paletteType)
    {
    case Palette_Jet:
        {
            if (x < 0.0) x = 0.0;
            else if (x > 1.0) x = 1.0;
            x *= numPalEntries;
            int n = (int) x;
            return paletteDataJet[n];
        }
        break;
    case Palette_Autumn:
        {
            if (x < 0.0) x = 0.0;
            else if (x > 1.0) x = 1.0;
            x *= numPalEntries;
            int n = (int) x;
            return paletteDataAutumn[n];
        }
        break;
    case Palette_Copper:
        {
            if (x < 0.0) x = 0.0;
            else if (x > 1.0) x = 1.0;
            x *= numPalEntries;
            int n = (int) x;
            return paletteDataCopper[n];
        }
        break;
    case Palette_Hot:
        {
            if (x < 0.0) x = 0.0;
            else if (x > 1.0) x = 1.0;
            x *= numPalEntries;
            int n = (int) x;
            return paletteDataHot[n];
        }
        break;
    case Palette_Cool:
        {
            if (x < 0.0) x = 0.0;
            else if (x > 1.0) x = 1.0;
            x *= numPalEntries;
            int n = (int) x;
            return paletteDataCool[n];
        }
        break;
    case Palette_BWAsc:
        {
            static float color[3];
            color[0] = color[1] = color[2] = x;
            return color;
        }
        break;
    case Palette_BWDesc:
        {
            static float color[3];
            color[0] = color[1] = color[2] = 1.0 - x;
            return color;
        }
        break;
    }
}

void SceneView::paletteCreate()
{
    int i;
    unsigned char palette[256][3];
    for (i = 0; i < Util::config()->paletteSteps; i++)
    {
        const float* color = paletteColor((double) i / Util::config()->paletteSteps);
        palette[i][0] = (unsigned char) (color[0] * 255);
        palette[i][1] = (unsigned char) (color[1] * 255);
        palette[i][2] = (unsigned char) (color[2] * 255);
    }
    for (i = Util::config()->paletteSteps; i < 256; i++)
        memcpy(palette[i], palette[Util::config()->paletteSteps-1], 3);

    glBindTexture(GL_TEXTURE_1D, 1);
    glTexImage1D(GL_TEXTURE_1D, 0, 3, 256, 0, GL_RGB, GL_UNSIGNED_BYTE, palette);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL);
}

void SceneView::paletteFilter()
{
    int palFilter = Util::config()->paletteFilter ? GL_LINEAR : GL_NEAREST;
    glBindTexture(GL_TEXTURE_1D, 1);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, palFilter);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, palFilter);
    paletteUpdateTexAdjust();
}

void SceneView::paletteUpdateTexAdjust()
{
    if (Util::config()->paletteFilter)
    {
        m_texScale = (double) (Util::config()->paletteSteps-1) / 256.0;
        m_texShift = 0.5 / 256.0;
    }
    else
    {
        m_texScale = (double) Util::config()->paletteSteps / 256.0;
        m_texShift = 0.0;
    }
}

void SceneView::initLighting()
{
    if (Util::config()->scalarView3DLighting || m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_Model)
    {
        glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_MODULATE);

        float light_specular[] = {  0.3f, 0.3f, 0.3f, 1.0f };
        float light_ambient[]  = {  0.1f, 0.1f, 0.1f, 1.0f };
        float light_diffuse[]  = {  0.8f, 0.8f, 0.8f, 0.9f };
        float light_position[] = {  0.0f, 10.0f, 0.0f, 1.0f };

        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
        glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
        // glLightfv(GL_LIGHT0, GL_POSITION, light_position);

        float material_specular[] = { 0.5f, 0.5f, 0.5f, 0.5f };
        float material_ambient[]  = { 0.5f, 0.5f, 0.5f, 1.0f };
        float material_diffuse[]  = { 0.8f, 0.8f, 0.8f, 1.0f };

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 8);
        glDisable(GL_COLOR_MATERIAL);

        glShadeModel(GL_SMOOTH);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        #if defined(GL_LIGHT_MODEL_COLOR_CONTROL) && defined(GL_SEPARATE_SPECULAR_COLOR)
            glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
        #endif

        /*
        float light_specular[] = {  0.1f, 0.1f, 0.1f, 1.0f };
        float light_ambient[]  = {  0.1f, 0.1f, 0.1f, 1.0f };
        float light_diffuse[]  = {  0.8f, 0.8f, 0.8f, 0.9f };
        float light_position[] = {  -100.0f, -60.0f, 10.0f, 0.0f };

        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
        glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
        glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
        glLightfv(GL_LIGHT0, GL_POSITION, light_position);

        float material_ambient[]  = { 0.5f, 0.5f, 0.5f, 1.0f };
        float material_diffuse[]  = { 0.8f, 0.8f, 0.8f, 0.8f };
        float material_specular[] = { 0.5f, 0.5f, 0.5f, 0.5f };

        glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, material_ambient);
        glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, material_diffuse);
        glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, material_specular);
        glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 8);
        glDisable(GL_COLOR_MATERIAL);

        glShadeModel(GL_SMOOTH);
        glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, 1);
        #if defined(GL_LIGHT_MODEL_COLOR_CONTROL) && defined(GL_SEPARATE_SPECULAR_COLOR)
            glLightModeli(GL_LIGHT_MODEL_COLOR_CONTROL, GL_SEPARATE_SPECULAR_COLOR);
        #endif
        */
    }
}

// events *****************************************************************************************************************************

void SceneView::keyPressEvent(QKeyEvent *event)
{
    if (is3DMode())
    {

    }
    else
    {
        Point stepTemp = position(Point(contextWidth(), contextHeight()));
        stepTemp.x = stepTemp.x - m_offset2d.x;
        stepTemp.y = stepTemp.y - m_offset2d.y;
        double step = qMin(stepTemp.x, stepTemp.y) / 10.0;

        switch (event->key())
        {
        case Qt::Key_Up:
            {
                m_offset2d.y += step;
                refresh();
            }
            break;
        case Qt::Key_Down:
            {
                m_offset2d.y -= step;
                refresh();
            }
            break;
        case Qt::Key_Left:
            {
                m_offset2d.x -= step;
                refresh();
            }
            break;
        case Qt::Key_Right:
            {
                m_offset2d.x += step;
                refresh();
            }
            break;
        case Qt::Key_Plus:
            {
                doZoomIn();
            }
            break;
        case Qt::Key_Minus:
            {
                doZoomOut();
            }
            break;
        case Qt::Key_Delete:
            {
                m_scene->deleteSelected();
            }
            break;
        case Qt::Key_Space:
            {
                doSceneObjectProperties();
            }
            break;
        case Qt::Key_Escape:
            {
                m_scene->selectNone();
                emit mousePressed();
                refresh();
            }
            break;
        default:
            QGLWidget::keyPressEvent(event);
        }

        // snap to grid
        m_snapToGrid = ((Util::config()->snapToGrid) && (event->modifiers() & Qt::ControlModifier) && (m_sceneMode == SceneMode_OperateOnNodes));

        // select all
        if ((event->modifiers() & Qt::ControlModifier) && (event->key() == Qt::Key_A))
        {
            if (m_sceneMode == SceneMode_Postprocessor)
            {
                // select volume integral area
                if (actPostprocessorModeVolumeIntegral->isChecked())
                {
                    m_scene->selectAll(SceneMode_OperateOnLabels);
                    emit mousePressed();
                }

                // select surface integral area
                if (actPostprocessorModeSurfaceIntegral->isChecked())
                {
                    m_scene->selectAll(SceneMode_OperateOnEdges);
                    emit mousePressed();
                }
            }
            else
            {
                m_scene->selectAll(m_sceneMode);
            }

            refresh();
        }

        // add node with coordinates under mouse pointer
        if ((event->modifiers() & Qt::AltModifier & Qt::ControlModifier) | (event->key() == Qt::Key_N))
        {
            Point p = position(Point(m_lastPos.x(), m_lastPos.y()));
            m_scene->doNewNode(p);
        }
        if ((event->modifiers() & Qt::AltModifier & Qt::ControlModifier) | (event->key() == Qt::Key_L))
        {
            Point p = position(Point(m_lastPos.x(), m_lastPos.y()));
            m_scene->doNewLabel(p);
        }
    }
}

void SceneView::keyReleaseEvent(QKeyEvent *event)
{
    if (is3DMode())
    {

    }
    else
    {
        if (m_snapToGrid)
        {
            m_snapToGrid = false;
            updateGL();
        }
    }
}

void SceneView::mousePressEvent(QMouseEvent *event)
{
    m_lastPos = event->pos();

    if (is3DMode())
    {

    }
    else
    {
        Point p = position(Point(event->pos().x(), event->pos().y()));

        if (event->button() & Qt::LeftButton)
        {
            // zoom region
            if (actSceneZoomRegion->isChecked())
            {
                m_regionPos = m_lastPos;
                actSceneZoomRegion->setChecked(false);
                actSceneZoomRegion->setData(true);
                m_region = true;

                return;
            }

            // select region
            if (actSceneViewSelectRegion->isChecked())
            {
                m_regionPos = m_lastPos;
                actSceneViewSelectRegion->setChecked(false);
                actSceneViewSelectRegion->setData(true);
                m_region = true;

                return;
            }

            if ((m_sceneMode == SceneMode_Postprocessor) &&
                !(m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_ScalarView3D ||
                  m_sceneViewSettings.postprocessorShow == SceneViewPostprocessorShow_ScalarView3DSolid))
            {
                // local point value
                if (actPostprocessorModeLocalPointValue->isChecked())
                    emit mousePressed(p);
                // select volume integral area
                if (actPostprocessorModeVolumeIntegral->isChecked())
                {
                    int index = m_scene->sceneSolution()->findTriangleInMesh(m_scene->sceneSolution()->meshInitial(), p);
                    if (index != -1)
                    {
                        //  find label marker
                        int labelIndex = m_scene->sceneSolution()->meshInitial()->get_element_fast(index)->marker;

                        m_scene->labels[labelIndex]->isSelected = !m_scene->labels[labelIndex]->isSelected;
                        updateGL();
                    }
                    emit mousePressed();
                }
                // select surface integral area
                if (actPostprocessorModeSurfaceIntegral->isChecked())
                {
                    //  find edge marker
                    SceneEdge *edge = findClosestEdge(p);

                    edge->isSelected = !edge->isSelected;
                    updateGL();

                    emit mousePressed();
                }
            }
        }

        // add node edge or label by mouse click
        if (event->modifiers() & Qt::ControlModifier)
        {
            // add node directly by mouse click
            if (m_sceneMode == SceneMode_OperateOnNodes)
            {
                Point pointNode;

                // snap to grid
                if (m_snapToGrid)
                {
                    Point snapPoint = position(Point(m_lastPos.x(), m_lastPos.y()));

                    pointNode.x = round(snapPoint.x / Util::config()->gridStep) * Util::config()->gridStep;
                    pointNode.y = round(snapPoint.y / Util::config()->gridStep) * Util::config()->gridStep;
                }
                else
                {
                    pointNode = p;
                }

                SceneNode *node = new SceneNode(pointNode);
                SceneNode *nodeAdded = m_scene->addNode(node);
                if (nodeAdded == node) m_scene->undoStack()->push(new SceneNodeCommandAdd(node->point));
                updateGL();
            }
            if (m_sceneMode == SceneMode_OperateOnEdges)
            {
                // add edge directly by mouse click
                SceneNode *node = findClosestNode(p);
                if (node)
                {
                    if (m_nodeLast == NULL)
                    {
                        m_nodeLast = node;
                        m_nodeLast->isSelected = true;
                    }
                    else
                    {
                        if (node != m_nodeLast)
                        {
                            SceneEdge *edge = new SceneEdge(m_nodeLast, node, m_scene->edgeMarkers[0], 0);
                            SceneEdge *edgeAdded = m_scene->addEdge(edge);
                            if (edgeAdded == edge) m_scene->undoStack()->push(new SceneEdgeCommandAdd(edge->nodeStart->point,
                                                                                                      edge->nodeEnd->point,
                                                                                                      edge->marker->name,
                                                                                                      edge->angle));
                        }

                        m_nodeLast->isSelected = false;
                        m_nodeLast = NULL;
                    }

                    updateGL();
                }
            }
            // add label directly by mouse click
            if (m_sceneMode == SceneMode_OperateOnLabels)
            {
                SceneLabel *label = new SceneLabel(p, m_scene->labelMarkers[0], 0, 0);
                SceneLabel *labelAdded = m_scene->addLabel(label);
                if (labelAdded == label) m_scene->undoStack()->push(new SceneLabelCommandAdd(label->point,
                                                                                             label->marker->name,
                                                                                             label->area,
                                                                                             label->polynomialOrder));
                updateGL();
            }
        }

        if ((event->modifiers() == 0) && (event->button() & Qt::LeftButton))
        {
            // select scene objects
            if (m_sceneMode == SceneMode_OperateOnNodes)
            {
                // select the closest node
                SceneNode *node = findClosestNode(p);
                if (node)
                {
                    node->isSelected = !node->isSelected;
                    updateGL();
                }
            }

            if (m_sceneMode == SceneMode_OperateOnEdges)
            {
                // select the closest label
                SceneEdge *edge = findClosestEdge(p);
                if (edge)
                {
                    edge->isSelected = !edge->isSelected;
                    updateGL();
                }
            }
            if (m_sceneMode == SceneMode_OperateOnLabels)
            {
                // select the closest label
                SceneLabel *label = findClosestLabel(p);
                if (label)
                {
                    label->isSelected = !label->isSelected;
                    updateGL();
                }
            }
        }
    }
}

void SceneView::mouseDoubleClickEvent(QMouseEvent * event)
{
    if (is3DMode())
    {

    }
    else
    {
        Point p = position(Point(event->pos().x(), event->pos().y()));

        // zoom best fit
        if (!(event->modifiers() & Qt::ControlModifier))
        {
            if ((event->buttons() & Qt::MidButton) || ((event->buttons() & Qt::LeftButton) && (event->modifiers() & Qt::ShiftModifier)))
            {
                doZoomBestFit();
            }

            if (event->button() & Qt::LeftButton)
            {
                // select scene objects
                m_scene->selectNone();
                if (m_sceneMode == SceneMode_OperateOnNodes)
                {
                    // select the closest node
                    SceneNode *node = findClosestNode(p);
                    if (node)
                    {
                        node->isSelected = true;
                        updateGL();
                        if (node->showDialog(this) == QDialog::Accepted)
                        {
                            updateGL();
                        }
                    }
                }
                if (m_sceneMode == SceneMode_OperateOnEdges)
                {
                    // select the closest label
                    SceneEdge *edge = findClosestEdge(p);
                    if (edge)
                    {
                        edge->isSelected = true;
                        updateGL();
                        if (edge->showDialog(this) == QDialog::Accepted)
                        {
                            updateGL();
                        }
                    }
                }
                if (m_sceneMode == SceneMode_OperateOnLabels)
                {
                    // select the closest label
                    SceneLabel *label = findClosestLabel(p);
                    if (label)
                    {
                        label->isSelected = true;
                        updateGL();
                        if (label->showDialog(this) == QDialog::Accepted)
                        {
                            updateGL();
                        }
                    }
                }
                m_scene->selectNone();
                updateGL();
            }
        }
    }
}

void SceneView::mouseReleaseEvent(QMouseEvent *event)
{
    setCursor(Qt::ArrowCursor);

    if (is3DMode())
    {

    }
    else
    {
        // zoom or select region
        actSceneZoomRegion->setChecked(false);
        actSceneViewSelectRegion->setChecked(false);

        if (m_region)
        {
            Point posStart = position(Point(m_regionPos.x(), m_regionPos.y()));
            Point posEnd = position(Point(m_lastPos.x(), m_lastPos.y()));

            if (actSceneZoomRegion->data().value<bool>())
                doZoomRegion(Point(qMin(posStart.x, posEnd.x), qMin(posStart.y, posEnd.y)), Point(qMax(posStart.x, posEnd.x), qMax(posStart.y, posEnd.y)));
            if (actSceneViewSelectRegion->data().value<bool>())
                selectRegion(Point(qMin(posStart.x, posEnd.x), qMin(posStart.y, posEnd.y)), Point(qMax(posStart.x, posEnd.x), qMax(posStart.y, posEnd.y)));

            actSceneZoomRegion->setData(false);
            actSceneViewSelectRegion->setData(false);
        }

        m_region = false;
    }

    updateGL();
}

void SceneView::mouseMoveEvent(QMouseEvent *event)
{
    int dx = event->x() - m_lastPos.x();
    int dy = event->y() - m_lastPos.y();

    m_lastPos = event->pos();

    setToolTip("");

    if (is3DMode())
    {
        // pan
        if ((event->buttons() & Qt::MidButton) || ((event->buttons() & Qt::LeftButton) && (event->modifiers() & Qt::ShiftModifier)))
        {
            setCursor(Qt::PointingHandCursor);

            m_offset3d.x -= 2.0/contextWidth() * dx*aspect();
            m_offset3d.y += 2.0/contextHeight() * dy;

            updateGL();
        }

        // rotate
        if ((event->buttons() & Qt::LeftButton) && !(event->modifiers() & Qt::ShiftModifier) && !(event->modifiers() & Qt::ControlModifier))
        {
            setCursor(Qt::PointingHandCursor);

            m_rotation3d.x -= dy;
            m_rotation3d.y += dx;

            updateGL();
        }
        if ((event->buttons() & Qt::LeftButton) && (event->modifiers() & Qt::ControlModifier))
        {
            setCursor(Qt::PointingHandCursor);

            m_rotation3d.z -= dy;

            updateGL();
        }
    }
    else
    {
        Point p = position(Point(m_lastPos.x(), m_lastPos.y()));

        // zoom or select region
        if (m_region)
            updateGL();

        // snap to grid
        if (m_snapToGrid && !(event->modifiers() & Qt::ControlModifier))
        {
            m_snapToGrid = false;
            updateGL();
        }
        m_snapToGrid = ((Util::config()->snapToGrid) && (event->modifiers() & Qt::ControlModifier) && (m_sceneMode == SceneMode_OperateOnNodes));

        if (m_snapToGrid)
            updateGL();

        // pan
        if ((event->buttons() & Qt::MidButton) || ((event->buttons() & Qt::LeftButton) && (event->modifiers() & Qt::ShiftModifier)))
        {
            setCursor(Qt::PointingHandCursor);

            m_offset2d.x -= 2.0/contextWidth() * dx/m_scale2d*aspect();
            m_offset2d.y += 2.0/contextHeight() * dy/m_scale2d;

            updateGL();
        }

        // hints
        if (event->modifiers() == 0)
        {
            // highlight scene objects
            if (m_sceneMode == SceneMode_OperateOnNodes)
            {
                // highlight the closest node
                SceneNode *node = findClosestNode(p);
                if (node)
                {
                    m_scene->highlightNone();
                    node->isHighlighted = true;
                    setToolTip(tr("<h3>Node</h3>Point: [%1; %2]<br/>Index: %3").
                               arg(node->point.x, 0, 'g', 3).
                               arg(node->point.y, 0, 'g', 3).
                               arg(m_scene->nodes.indexOf(node)));
                    updateGL();
                }
            }
            if (m_sceneMode == SceneMode_OperateOnEdges)
            {
                // highlight the closest label
                SceneEdge *edge = findClosestEdge(p);
                if (edge)
                {
                    m_scene->highlightNone();
                    edge->isHighlighted = true;
                    setToolTip(tr("<h3>Edge</h3>Point: [%1; %2] - [%3; %4]<br/>Boundary Condition: %5<br/>Angle: %6 deg.<br/>Index: %7 %8").
                               arg(edge->nodeStart->point.x, 0, 'g', 3).
                               arg(edge->nodeStart->point.y, 0, 'g', 3).
                               arg(edge->nodeEnd->point.x, 0, 'g', 3).
                               arg(edge->nodeEnd->point.y, 0, 'g', 3).
                               arg(edge->marker->name).
                               arg(edge->angle, 0, 'f', 0).
                               arg(m_scene->edges.indexOf(edge)).
                               arg(edge->marker->html()));
                    updateGL();
                }
            }
            if (m_sceneMode == SceneMode_OperateOnLabels)
            {
                // highlight the closest label
                SceneLabel *label = findClosestLabel(p);
                if (label)
                {
                    m_scene->highlightNone();
                    label->isHighlighted = true;
                    setToolTip(tr("<h3>Label</h3>Point: [%1; %2]<br/>Material: %3<br/>Triangle area: %4 m<sup>2</sup><br/>Polynomial order: %5<br/>Index: %6 %7").
                               arg(label->point.x, 0, 'g', 3).
                               arg(label->point.y, 0, 'g', 3).
                               arg(label->marker->name).
                               arg(label->area, 0, 'g', 3).
                               arg(label->polynomialOrder).
                               arg(m_scene->labels.indexOf(label)).
                               arg(label->marker->html()));
                    updateGL();
                }
            }
        }

        if (event->modifiers() & Qt::ControlModifier)
        {
            // add edge directly by mouse click - highlight
            if (m_sceneMode == SceneMode_OperateOnEdges)
            {
                // add edge directly by mouse click
                SceneNode *node = findClosestNode(p);
                if (node)
                {
                    m_scene->highlightNone();
                    node->isHighlighted = true;
                    updateGL();
                }
            }
        }


        if (m_snapToGrid)
        {
            Point snapPoint;
            snapPoint.x = round(p.x / Util::config()->gridStep) * Util::config()->gridStep;
            snapPoint.y = round(p.y / Util::config()->gridStep) * Util::config()->gridStep;

            emit mouseMoved(QPointF(snapPoint.x, snapPoint.y));
        }
        else
        {
            emit mouseMoved(QPointF(p.x, p.y));
        }
    }
}

void SceneView::wheelEvent(QWheelEvent *event)
{
    if (is3DMode())
    {
        setZoom(event->delta()/150.0);
    }
    else
    {
        if (Util::config()->zoomToMouse)
        {
            Point posMouse;
            posMouse = Point((2.0/contextWidth()*(event->pos().x() - contextWidth()/2.0))/m_scale2d*aspect(),
                            -(2.0/contextHeight()*(event->pos().y() - contextHeight()/2.0))/m_scale2d);

            m_offset2d.x += posMouse.x;
            m_offset2d.y += posMouse.y;

            m_scale2d = m_scale2d * pow(1.2, event->delta()/150.0);

            posMouse = Point((2.0/contextWidth()*(event->pos().x() - contextWidth()/2.0))/m_scale2d*aspect(),
                            -(2.0/contextHeight()*(event->pos().y() - contextHeight()/2.0))/m_scale2d);

            m_offset2d.x -= posMouse.x;
            m_offset2d.y -= posMouse.y;

            updateGL();
        }
        else
        {
            setZoom(event->delta()/150.0);
        }
    }
}

void SceneView::contextMenuEvent(QContextMenuEvent *event)
{
    actSceneObjectProperties->setEnabled(false);

    // set boundary context menu
    if (m_sceneMode == SceneMode_OperateOnEdges)
        actSceneObjectProperties->setEnabled(m_scene->selectedCount() > 0);

    // set material context menu
    if (m_sceneMode == SceneMode_OperateOnLabels)
        actSceneObjectProperties->setEnabled(m_scene->selectedCount() > 0);


    mnuInfo->exec(event->globalPos());
}

void SceneView::closeEvent(QCloseEvent *event)
{
    event->ignore();
}

// slots *****************************************************************************************************************************

void SceneView::doZoomBestFit()
{
    RectPoint rect = m_scene->boundingBox();
    doZoomRegion(rect.start, rect.end);
}

void SceneView::doZoomIn()
{
    setZoom(1.2);
}

void SceneView::doZoomOut()
{
    setZoom(-1/1.2);
}

void SceneView::doZoomRegion(const Point &start, const Point &end)
{
    if (fabs(end.x-start.x) < EPS_ZERO || fabs(end.y-start.y) < EPS_ZERO) return;

    m_offset2d.x = (start.x+end.x)/2.0;
    m_offset2d.y = (start.y+end.y)/2.0;

    double sceneWidth = end.x-start.x;
    double sceneHeight = end.y-start.y;

    double maxScene = (((double) contextWidth() / (double) contextHeight()) < (sceneWidth / sceneHeight)) ? sceneWidth/aspect() : sceneHeight;

    if (maxScene > 0.0)
        m_scale2d = 1.95/maxScene;

    setZoom(0);
}

void SceneView::doSetChartLine(const Point &start, const Point &end)
{
    // set line for chart
    m_chartLine.start = start;
    m_chartLine.end = end;

    updateGL();
}

void SceneView::doDefaultValues()
{
    m_snapToGrid = false;
    m_region = false;
    m_isSolutionPrepared = false;

    // 2d
    m_scale2d = 1.0;
    m_offset2d = Point();

    // 3d
    m_scale3d = 1.0;
    m_offset3d = Point();
    m_rotation3d = Point3();

    m_chartLine.start = Point();
    m_chartLine.end = Point();

    m_sceneViewSettings.defaultValues();

    doInvalidated();
    doZoomBestFit();

    actPostprocessorModeLocalPointValue->trigger();
}

void SceneView::solved()
{
    actSceneModePostprocessor->trigger();
    m_sceneViewSettings.showInitialMesh = false;
}

void SceneView::doInvalidated()
{
    resizeGL(width(), height());

    if (!m_scene->sceneSolution()->isSolved())
    {
        if (m_sceneMode == SceneMode_Postprocessor)
            actSceneModeNode->trigger();
    }

    actSceneModePostprocessor->setEnabled(m_scene->sceneSolution()->isSolved());
    actSceneViewSelectMarker->setEnabled(m_scene->sceneSolution()->isSolved());
    actSceneZoomRegion->setEnabled(!is3DMode());

    m_scene->actDeleteSelected->setEnabled(m_sceneMode != SceneMode_Postprocessor);
    actSceneViewSelectRegion->setEnabled(m_sceneMode != SceneMode_Postprocessor);

    actPostprocessorModeLocalPointValue->setEnabled(m_sceneMode == SceneMode_Postprocessor && !is3DMode());
    actPostprocessorModeSurfaceIntegral->setEnabled(m_sceneMode == SceneMode_Postprocessor && !is3DMode());
    actPostprocessorModeVolumeIntegral->setEnabled(m_sceneMode == SceneMode_Postprocessor && !is3DMode());

    emit mousePressed();

    paintGL();
    updateGL();
}

void SceneView::timeStepChanged(bool showViewProgress)
{
    qDebug() << "SceneView::timeStepChanged()";

    if (!Util::scene()->sceneSolution()->isSolving())
    {
        if (showViewProgress)
        {
            ProgressDialog progressDialog;
            progressDialog.appendProgressItem(new ProgressItemProcessView());
            progressDialog.run();
        }
        else
        {
            ProgressItemProcessView progressItemProcessView;
            progressItemProcessView.run();
        }
    }

    clearGLLists();
    m_isSolutionPrepared = true;

    refresh();
}

void SceneView::refresh()
{
    paintGL();
    updateGL();
}

void SceneView::doMaterialGroup(QAction *action)
{
    if (SceneLabelMarker *labelMarker = action->data().value<SceneLabelMarker *>())
        m_scene->setLabelLabelMarker(labelMarker);
}

void SceneView::doBoundaryGroup(QAction *action)
{
    if (SceneEdgeMarker *edgeMarker = action->data().value<SceneEdgeMarker *>())
        m_scene->setEdgeEdgeMarker(edgeMarker);
}

void SceneView::doShowGroup(QAction *action)
{
    m_sceneViewSettings.showContours = actShowContours->isChecked();
    m_sceneViewSettings.showVectors = actShowVectors->isChecked();
    m_sceneViewSettings.showSolutionMesh = actShowSolutionMesh->isChecked();

    doInvalidated();
}

void SceneView::doPostprocessorModeGroup(QAction *action)
{
    m_scene->selectNone();
    updateGL();
}

void SceneView::doSceneViewProperties()
{
    SceneViewPostprocessorShow postprocessorShow = m_sceneViewSettings.postprocessorShow;

    SceneViewDialog sceneViewDialog(this, this);
    if (sceneViewDialog.showDialog() == QDialog::Accepted)
    {
        // set defaults
        if (postprocessorShow != m_sceneViewSettings.postprocessorShow)
        {
            if (is3DMode())
            {
                m_rotation3d.x = 66.0;
                m_rotation3d.y = -35.0;
                m_rotation3d.z = 0.0;

                m_offset3d.x = 0.0;
                m_offset3d.y = 0.0;

                m_scale3d = 0.6 * m_scale2d;
            }

            doZoomBestFit();
        }
    }
}

void SceneView::doSceneObjectProperties()
{
    if (m_sceneMode == SceneMode_OperateOnEdges)
    {
        if (m_scene->selectedCount() > 1)
        {
            EdgeMarkerDialog edgeMarkerDialog(this);
            edgeMarkerDialog.exec();
        }
        if (m_scene->selectedCount() == 1)
        {
            for (int i = 0; i < m_scene->edges.count(); i++)
            {
                if (m_scene->edges[i]->isSelected)
                    m_scene->edges[i]->showDialog(this);
            }
        }
    }
    if (m_sceneMode == SceneMode_OperateOnLabels)
    {
        if (m_scene->selectedCount() > 1)
        {
            LabelMarkerDialog labelMarkerDialog(this);
            labelMarkerDialog.exec();
        }
        if (m_scene->selectedCount() == 1)
        {
            for (int i = 0; i < m_scene->labels.count(); i++)
            {
                if (m_scene->labels[i]->isSelected)
                    m_scene->labels[i]->showDialog(this);
            }
        }
    }

    m_scene->selectNone();
}

void SceneView::doSceneModeSet(QAction *action)
{
    actSceneModePostprocessor->setEnabled(m_scene->sceneSolution()->isSolved());

    if (actSceneModeNode->isChecked()) m_sceneMode = SceneMode_OperateOnNodes;
    if (actSceneModeEdge->isChecked()) m_sceneMode = SceneMode_OperateOnEdges;
    if (actSceneModeLabel->isChecked()) m_sceneMode = SceneMode_OperateOnLabels;
    if (actSceneModePostprocessor->isChecked()) m_sceneMode = SceneMode_Postprocessor;

    m_scene->highlightNone();
    m_scene->selectNone();

    switch (m_sceneMode)
    {
    case SceneMode_OperateOnNodes:
        break;
    case SceneMode_OperateOnEdges:
        m_nodeLast = NULL;
        break;
    case SceneMode_OperateOnLabels:
        break;
    case SceneMode_Postprocessor:
        break;
    }

    doInvalidated();

    emit sceneModeChanged(m_sceneMode);
}

void SceneView::doSelectMarker()
{
    SceneMarkerSelectDialog sceneMarkerSelectDialog(this, QApplication::activeWindow());
    sceneMarkerSelectDialog.exec();
}

void SceneView::doSelectBasic()
{
    SceneBasicSelectDialog sceneBasicSelectDialog(this, QApplication::activeWindow());
    sceneBasicSelectDialog.exec();
}

void SceneView::processedRangeContour()
{

}

void SceneView::processedRangeScalar()
{
    paletteFilter();
    paletteUpdateTexAdjust();
    paletteCreate();

    if (m_sceneViewSettings.scalarRangeAuto)
    {
        m_sceneViewSettings.scalarRangeMin = m_scene->sceneSolution()->linScalarView().get_min_value();
        m_sceneViewSettings.scalarRangeMax = m_scene->sceneSolution()->linScalarView().get_max_value();
    }
}

void SceneView::processedRangeVector()
{
}

void SceneView::setZoom(double power)
{
    if (is3DMode())
    {
        m_scale3d = m_scale3d * pow(1.2, power);
    }
    else
    {
        m_scale2d = m_scale2d * pow(1.2, power);
    }

    updateGL();

    Point p(pos().x(), pos().y());
    emit mouseMoved(QPointF(position(p).x, position(p).y));
}

void SceneView::selectRegion(const Point &start, const Point &end)
{
    m_scene->selectNone();

    switch (m_sceneMode)
    {
    case SceneMode_OperateOnNodes:
        foreach (SceneNode *node, m_scene->nodes)
            if (node->point.x >= start.x && node->point.x <= end.x && node->point.y >= start.y && node->point.y <= end.y)
                node->isSelected = true;
        break;
    case SceneMode_OperateOnEdges:
        foreach (SceneEdge *edge, m_scene->edges)
            if (edge->nodeStart->point.x >= start.x && edge->nodeStart->point.x <= end.x && edge->nodeStart->point.y >= start.y && edge->nodeStart->point.y <= end.y &&
                edge->nodeEnd->point.x >= start.x && edge->nodeEnd->point.x <= end.x && edge->nodeEnd->point.y >= start.y && edge->nodeEnd->point.y <= end.y)
                edge->isSelected = true;
        break;
    case SceneMode_OperateOnLabels:
        foreach (SceneLabel *label, m_scene->labels)
            if (label->point.x >= start.x && label->point.x <= end.x && label->point.y >= start.y && label->point.y <= end.y)
                label->isSelected = true;
        break;
    }
}

SceneNode *SceneView::findClosestNode(const Point &point)
{
    SceneNode *nodeClosest = NULL;

    double distance = CONST_DOUBLE;
    foreach (SceneNode *node, m_scene->nodes)
    {
        double nodeDistance = node->distance(point);
        if (node->distance(point) < distance)
        {
            distance = nodeDistance;
            nodeClosest = node;
        }
    }

    return nodeClosest;
}

SceneEdge *SceneView::findClosestEdge(const Point &point)
{
    SceneEdge *edgeClosest = NULL;

    double distance = CONST_DOUBLE;
    foreach (SceneEdge *edge, m_scene->edges)
    {
        double edgeDistance = edge->distance(point);
        if (edge->distance(point) < distance)
        {
            distance = edgeDistance;
            edgeClosest = edge;
        }
    }

    return edgeClosest;
}

SceneLabel *SceneView::findClosestLabel(const Point &point)
{
    SceneLabel *labelClosest = NULL;

    double distance = CONST_DOUBLE;
    foreach (SceneLabel *label, m_scene->labels)
    {
        double labelDistance = label->distance(point);
        if (label->distance(point) < distance)
        {
            distance = labelDistance;
            labelClosest = label;
        }
    }

    return labelClosest;
}

void SceneView::drawArc(const Point &point, double r, double startAngle, double arcAngle, int segments)
{
    double theta = arcAngle / double(segments - 1);

    glBegin(GL_LINE_STRIP);
    for (int i = 0; i < segments; i++)
    {
        double arc = (startAngle + i*theta)/180.0*M_PI;

        double x = r * cos(arc);
        double y = r * sin(arc);

        glVertex3d(point.x + x, point.y + y, 0.0);
    }
    glEnd();
}

void SceneView::drawBlend(Point start, Point end, double red, double green, double blue, double alpha)
{
    // store color
    double color[4];
    glGetDoublev(GL_CURRENT_COLOR, color);

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    // blended rectangle
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4f(red, green, blue, alpha);

    glBegin(GL_QUADS);
    glVertex2d(start.x, start.y);
    glVertex2d(end.x, start.y);
    glVertex2d(end.x, end.y);
    glVertex2d(start.x, end.y);
    glEnd();

    glDisable(GL_BLEND);
    glDisable(GL_POLYGON_OFFSET_FILL);

    // retrieve color
    glColor4d(color[0], color[1], color[2], color[3]);
}

void SceneView::renderTextPos(double x, double y, double z, const QString &str, bool blend)
{
    if (blend)
    {
        Point size((2.0/contextWidth()*fontMetrics().width(" "))/m_scale2d*aspect(),
                   (2.0/contextHeight()*fontMetrics().height())/m_scale2d);

        double xs = x - size.x / 2.0;
        double ys = y - size.y * 1.15 / 3.2;
        double xe = xs + size.x * (str.size() + 1);
        double ye = ys + size.y * 1.15;

        drawBlend(Point(xs, ys), Point(xe, ye));
    }

    renderText(x, y, z, str);
}

void SceneView::paintPostprocessorSelectedVolume()
{
    if (!m_scene->sceneSolution()->isMeshed()) return;

    glEnable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glColor4d(Util::config()->colorSelected.redF(),
              Util::config()->colorSelected.greenF(),
              Util::config()->colorSelected.blueF(),
              0.5);

    // triangles
    glBegin(GL_TRIANGLES);
    for (int i = 0; i < m_scene->sceneSolution()->meshInitial()->get_num_active_elements(); i++)
    {
        Element *element = m_scene->sceneSolution()->meshInitial()->get_element(i);
        if (m_scene->labels[element->marker]->isSelected)
        {
            if (element->is_triangle())
            {
                glVertex2d(element->vn[0]->x, element->vn[0]->y);
                glVertex2d(element->vn[1]->x, element->vn[1]->y);
                glVertex2d(element->vn[2]->x, element->vn[2]->y);
            }
        }
    }
    glEnd();

    glDisable(GL_BLEND);
    glDisable(GL_POLYGON_OFFSET_FILL);
}

void SceneView::paintPostprocessorSelectedSurface()
{
    // edges
    foreach (SceneEdge *edge, m_scene->edges) {
        glColor3d(Util::config()->colorSelected.redF(), Util::config()->colorSelected.greenF(), Util::config()->colorSelected.blueF());
        glLineWidth(3.0);

        if (edge->isSelected)
        {
            if (edge->angle == 0)
            {
                glBegin(GL_LINES);
                glVertex2d(edge->nodeStart->point.x, edge->nodeStart->point.y);
                glVertex2d(edge->nodeEnd->point.x, edge->nodeEnd->point.y);
                glEnd();
            }
            else
            {
                Point center = edge->center();
                double radius = edge->radius();
                double startAngle = atan2(center.y - edge->nodeStart->point.y, center.x - edge->nodeStart->point.x) / M_PI*180 - 180;

                drawArc(center, radius, startAngle, edge->angle, edge->angle/5);
            }
        }
        glLineWidth(1.0);
    }
}

ErrorResult SceneView::saveImageToFile(const QString &fileName, int w, int h)
{
    QPixmap pixmap = renderScenePixmap(w, h);
    if (pixmap.save(fileName, "PNG"))
        resizeGL(width(), height());
    else
        return ErrorResult(ErrorResultType_Critical, tr("Image cannot be saved to the file '%1'.").arg(fileName));

    return ErrorResult();
}

void SceneView::saveImagesForReport(const QString &path, bool showRulers, bool showGrid, int w, int h)
{
    // store sceneview settings
    SceneViewSettings sceneViewSettingsCopy = m_sceneViewSettings;
    SceneMode sceneModeCopy = m_sceneMode;
    double scale2dCopy = m_scale2d;
    Point offset2dCopy = m_offset2d;
    Point offset3dCopy = m_offset3d;
    Point3 rotation3dCopy = m_rotation3d;

    bool showRulersCopy = Util::config()->showRulers;
    bool showGridCopy = m_sceneViewSettings.showGrid;

    // remove old files
    QFile::remove(path + "/geometry.png");
    QFile::remove(path + "/mesh.png");
    QFile::remove(path + "/order.png");
    QFile::remove(path + "/scalarview.png");

    doZoomBestFit();

    m_sceneViewSettings.showGeometry = true;
    m_sceneViewSettings.showContours = false;
    m_sceneViewSettings.showVectors = false;
    m_sceneViewSettings.showInitialMesh = false;
    m_sceneViewSettings.showSolutionMesh = false;

    Util::config()->showRulers = showRulers;
    m_sceneViewSettings.showGrid = showGrid;

    // geometry
    actSceneModeLabel->trigger();
    ErrorResult resultGeometry = saveImageToFile(path + "/geometry.png", w, h);
    if (resultGeometry.isError())
        resultGeometry.showDialog();

    // mesh
    if (m_scene->sceneSolution()->isMeshed())
    {
        // show only initial mesh
        actSceneModeLabel->trigger();

        m_sceneViewSettings.showInitialMesh = true;
        ErrorResult resultMesh1 = saveImageToFile(path + "/mesh.png", w, h);
        if (resultMesh1.isError())
            resultMesh1.showDialog();
        m_sceneViewSettings.showInitialMesh = false;
    }
    if (m_scene->sceneSolution()->isSolved())
    {
        // when solved show both meshes
        actSceneModePostprocessor->trigger();

        m_sceneViewSettings.postprocessorShow = SceneViewPostprocessorShow_None;
        updateGL();

        m_sceneViewSettings.showInitialMesh = true;
        m_sceneViewSettings.showSolutionMesh = true;
        ErrorResult resultMesh2 = saveImageToFile(path + "/mesh.png", w, h);
        if (resultMesh2.isError())
            resultMesh2.showDialog();
        m_sceneViewSettings.showInitialMesh = false;
        m_sceneViewSettings.showSolutionMesh = false;
    }

    // order
    m_sceneViewSettings.postprocessorShow = SceneViewPostprocessorShow_Order;
    updateGL();
    ErrorResult resultOrder = saveImageToFile(path + "/order.png", w, h);
    if (resultOrder.isError())
        resultOrder.showDialog();

    if (m_scene->sceneSolution()->isSolved())
    {
        actSceneModePostprocessor->trigger();

        // last step
        if (m_scene->problemInfo()->hermes()->hasTransient())
            m_scene->sceneSolution()->setTimeStep(m_scene->sceneSolution()->timeStepCount() - 1);

        m_sceneViewSettings.scalarRangeAuto = true;
        m_sceneViewSettings.scalarPhysicFieldVariable = m_scene->problemInfo()->hermes()->scalarPhysicFieldVariable();
        m_sceneViewSettings.scalarPhysicFieldVariableComp = m_scene->problemInfo()->hermes()->scalarPhysicFieldVariableComp();
        m_sceneViewSettings.vectorPhysicFieldVariable = m_scene->problemInfo()->hermes()->vectorPhysicFieldVariable();

        doInvalidated();

        // scalar field
        m_sceneViewSettings.postprocessorShow = SceneViewPostprocessorShow_ScalarView;
        updateGL();
        ErrorResult resultScalarView = saveImageToFile(path + "/scalarview.png", w, h);
        if (resultScalarView.isError())
            resultScalarView.showDialog();
    }

    // restore sceneview settings
    m_sceneViewSettings = sceneViewSettingsCopy;
    m_sceneMode = sceneModeCopy;
    m_scale2d = scale2dCopy;
    m_offset2d = offset2dCopy;
    m_offset3d = offset3dCopy;
    m_rotation3d = rotation3dCopy;

    Util::config()->showRulers = showRulersCopy;
    m_sceneViewSettings.showGrid = showGridCopy;

    if (m_sceneMode == SceneMode_OperateOnNodes) actSceneModeNode->trigger();
    if (m_sceneMode == SceneMode_OperateOnLabels) actSceneModeEdge->isChecked();
    if (m_sceneMode == SceneMode_OperateOnLabels) actSceneModeLabel->isChecked();
    if (m_sceneMode == SceneMode_Postprocessor) actSceneModePostprocessor->isChecked();

    refresh();
}
