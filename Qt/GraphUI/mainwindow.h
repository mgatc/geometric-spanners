#ifndef MAINWINDOW_H
#define MAINWINDOW_H

#include <QMainWindow>
#include <QGraphicsScene>
#include <QGraphicsView>

#include "graphview.h"
#include "graphscene.h"
#include "playbackparameters.h"

QT_BEGIN_NAMESPACE
namespace Ui { class MainWindow; }
QT_END_NAMESPACE

class MainWindow : public QMainWindow
{
    Q_OBJECT

public:
    MainWindow(QWidget *parent = nullptr);
    ~MainWindow() override;

    void showEvent( QShowEvent* event ) override;

    GraphScene *_scene;
    GraphView *_view;
    PlaybackParameters *_controls;

private slots:
    void newFile();
    void open();
    void save();
    void exportPdf();
    void chooseAlgorithm();
    void help();
    void about();

private:
    Ui::MainWindow *ui;

    void createActions();
    void createMenus();

    QMenu *fileMenu;
    QMenu *algorithmMenu;
    QMenu *helpMenu;

    QAction *newAct;
    QAction *openAct;
    QAction *saveAct;
    QAction *exportPdfAct;
    QAction *chooseAlgorithmAct;
    QAction *helpAct;
    QAction *aboutAct;

    QColor _VIEWING_PLANE_BG_COLOR = Qt::white;
};
#endif // MAINWINDOW_H
