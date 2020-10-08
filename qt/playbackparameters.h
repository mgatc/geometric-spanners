#ifndef PLAYBACKPARAMETERS_H
#define PLAYBACKPARAMETERS_H

#include <QDialog>

QT_BEGIN_NAMESPACE
class QAbstractButton;
class QCheckBox;
class QDialogButtonBox;
class QLabel;
class QSlider;
class QLineEdit;
class QTableWidget;
class QTextEdit;
class QWidget;
QT_END_NAMESPACE


class PlaybackParameters : public QDialog
{
    Q_OBJECT
public:
    explicit PlaybackParameters(QWidget *parent = nullptr);

    bool playing;
    bool pausing;
    enum { SpanningGraph, PolygonSpanner, PlanarSpanner } algorithm;
    int framesPerSecond;
    long position;
    long length;
    int zoom;

public slots:
    void rcv_playRequest();           // from play button
    void rcv_startPauseRequest();
    void rcv_endPauseRequest();
    void rcv_playingInfo( bool playing ); // from GraphView
    void rcv_scrubRequest( int value ); // from slider
    void rcv_lengthInfo( int length ); // from GraphView
    void rcv_positionInfo( int position ); // from GraphView
    //void rcv_fpsRequestInfo( int fps ); // from dropdown
//    void rcv_fpsInfo( int fps ); // from GraphView
//    void rcv_zoomRequestInfo( bool zoomIn );
    void rcv_stepBackRequest();
    void rcv_stepForwardRequest();

signals:
    void send_playingIcon( QIcon icon ); // connect to setIcon slot of playButton
    void send_positionInfo( int position ); // connect to setValue slot of positionSlider
    void send_lengthInfo( int position ); // connect to setMaximum slot of positionSlider
    void send_fpsInfo( int fps ); // connect to setValue slot of fps
    void send_zoomLevelInfo( double level ); // connect to setValue slot of zoom label
    void send_stepRequest( bool forward );
    void send_playRequest( bool play ); // connect to rcv_playRequest of GraphView
    void send_scrubRequest( int position ); // connect to rcv_scrubRequest of GraphView
    void send_zoomRequest( double level ); // connect to setZoom of GraphView
    void send_fpsRequest( int fps ); // connect to setFramerate of GraphView

private:
    QAbstractButton* playButton;
    QAbstractButton* stepBackButton;
    QAbstractButton* stepForwardButton;
    QAbstractButton* zoomOutButton;
    QAbstractButton* zoomInButton;
    QLabel* framesPerSecondLabel; // playback speed control
    QSlider* framesPerSecondSlider;
    QSlider* positionSlider;
    QDialogButtonBox* zoomButton;

};

#endif // PLAYBACKPARAMETERS_H
