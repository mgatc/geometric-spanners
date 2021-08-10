#include "mainwindow.h"
//#include "./ui_mainwindow.h"

#include <fstream>
#include <string>
#include <vector>

#include <QTimer>
#include <QtWidgets>

#include "graphscene.h"
#include "graphview.h"
#include "tools/DelaunayGraph.h"#include "PlanarSpanner.h"
#include "playbackparameters.h"


MainWindow::MainWindow(QWidget *parent)
    : QMainWindow(parent)
    //, ui(new Ui::MainWindow)
{
    //ui->setupUi(this);

    QWidget *widget = new QWidget;
    setCentralWidget( widget );

    _scene = new GraphScene(this);

    _view = new GraphView( _scene, this );
    _view->scale(1, -1);
    _view->setBackgroundBrush( _VIEWING_PLANE_BG_COLOR );

    _controls = new PlaybackParameters(this);

    QVBoxLayout *layout = new QVBoxLayout;
    layout->setContentsMargins( 5, 5, 5, 5 );
    layout->addWidget( _view );
    layout->addWidget( _controls );

    widget->setLayout(layout);

    createActions();
    createMenus();

    setWindowTitle(tr("GraphUI"));
    setMinimumSize(600, 600);
    resize(960, 960);
    open();
}

MainWindow::~MainWindow()
{
    //delete ui;
}

void MainWindow::showEvent( QShowEvent* event )
{
    _view->fitInView( 0, 0, _scene->width(), _scene->height() );
    QWidget::showEvent(event);
}

void MainWindow::createActions()
{
    newAct = new QAction( tr( "&New" ) );
    newAct->setShortcuts( QKeySequence::New );
    newAct->setStatusTip( tr( "Create a new pointset" ) );
    connect( newAct, &QAction::triggered, this, &MainWindow::newFile );

    openAct = new QAction( tr( "&Open" ) );
    openAct->setShortcuts( QKeySequence::Open );
    openAct->setStatusTip( tr( "Open an existing pointset" ) );
    connect( openAct, &QAction::triggered, this, &MainWindow::open );

    saveAct = new QAction( tr( "&Save" ) );
    saveAct->setShortcuts( QKeySequence::Save );
    saveAct->setStatusTip( tr( "Save the open pointset" ) );
    connect( saveAct, &QAction::triggered, this, &MainWindow::save );

    exportPdfAct = new QAction( tr( "&Export PDF" ) );
    exportPdfAct->setStatusTip( tr( "Export as PDF" ) );
    connect( exportPdfAct, &QAction::triggered, this, &MainWindow::exportPdf );

    chooseAlgorithmAct = new QAction( tr( "&Choose Algorithm" ) );
    chooseAlgorithmAct->setStatusTip( tr( "Choose the algorithm to perform" ) );
    connect( chooseAlgorithmAct, &QAction::triggered, this, &MainWindow::chooseAlgorithm );

    helpAct = new QAction( tr( "&Help" ) );
    helpAct->setStatusTip( tr( "Get help" ) );
    connect( helpAct, &QAction::triggered, this, &MainWindow::help );

    aboutAct = new QAction( tr( "&About" ) );
    aboutAct->setStatusTip( tr( "About..." ) );
    connect( aboutAct, &QAction::triggered, this, &MainWindow::about );

    // connect the play button to the view, and setup feedback loop
    connect( _controls, SIGNAL(send_playRequest(bool)), _view, SLOT(rcv_playRequest(bool)) );
    connect( _view, SIGNAL(send_playingInfo(bool)), _controls, SLOT(rcv_playingInfo(bool)) );
    // connect position slider to the view for scrubbing and playback position feedback
    connect( _controls, SIGNAL(send_scrubRequest(int)), _view, SLOT(setPosition(int)) );
    connect( _view, SIGNAL(send_positionInfo(int)), _controls, SLOT(rcv_positionInfo(int)) );
    // connect step buttons
    connect( _controls, SIGNAL(send_stepRequest(bool)), _view, SLOT(rcv_stepRequest(bool)) );


    // Don't connect anything to the _scene here, as the scene is destroyed when changing the point set
    // Put in open()
}

void MainWindow::createMenus()
{
    fileMenu = menuBar()->addMenu( tr( "&File" ) );
    fileMenu->addAction( newAct );
    fileMenu->addAction( openAct );
    fileMenu->addSeparator();
    fileMenu->addAction( saveAct );
    fileMenu->addSeparator();
    fileMenu->addAction( exportPdfAct );

    algorithmMenu = menuBar()->addMenu( tr( "&Algorithm" ) );
    algorithmMenu->addAction( chooseAlgorithmAct );

    helpMenu = menuBar()->addMenu( tr( "&Help" ) );
    helpMenu->addAction( helpAct );
    helpMenu->addAction( aboutAct );
}

void MainWindow::newFile()
{
    // prompt to save current file if unsaved
    // reset graph object
    // display blank canvas
}

void MainWindow::open()
{
    // prompt to save current file if unsaved

    // open file selection menu
    QFileDialog dialog(this);

    QStringList fileNames;

    // parse points
    if( dialog.exec() ) {
        delete _scene;

        _scene = new GraphScene(this);
        _scene->setPointset( dialog.selectedFiles() );
        _scene->setSceneRect( _scene->itemsBoundingRect() );

        _view->setScene( _scene );

        connect( _scene, SIGNAL(send_playRequest(bool)), _view, SLOT(rcv_playRequest(bool)) );
        connect( _view, SIGNAL(signalAdvance(int)), _scene, SLOT(advance(int)) );

        // setup communication of the overall number of frames
        connect( _scene, SIGNAL(send_lengthInfo(int)), _view, SLOT(rcv_lengthInfo(int)) );
        connect( _view, SIGNAL(send_lengthInfo(int)), _controls, SLOT(rcv_lengthInfo(int)) );

        _view->fitInView( _scene->sceneRect(), Qt::KeepAspectRatio ); //setScale(5);
        _view->scaleItemSize(0.04);

    }
}


void MainWindow::save()
{
    // save edge list
}

void MainWindow::exportPdf()
{
    // use GeometricSpannerPrinter
}

// VITAL
void MainWindow::chooseAlgorithm()
{
    // display algorithm selection window (include algorithm params)
    // for now, just call algorithm directly from here
    _scene->runAlgorithm( GraphScene::GraphAlgorithm::PlanarSpanner );
}

void MainWindow::help()
{
    // display help window
}

void MainWindow::about()
{
    // display about window
}
