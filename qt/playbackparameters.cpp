#include <iostream>
#include <QtWidgets>
#include "playbackparameters.h"

PlaybackParameters::PlaybackParameters(QWidget *parent) : QDialog(parent) {

    positionSlider = new QSlider( Qt::Horizontal, this );
    positionSlider->setSingleStep(1);
    connect( positionSlider, SIGNAL(sliderMoved(int)), this, SLOT(rcv_scrubRequest(int)) );
    connect( positionSlider, SIGNAL(sliderPressed()), this, SLOT(rcv_startPauseRequest()) );
    connect( positionSlider, SIGNAL(sliderReleased()), this, SLOT(rcv_endPauseRequest()) );

    playButton = new QToolButton(this);
    playButton->setIcon( style()->standardIcon(QStyle::SP_MediaPlay) );
    connect( playButton, SIGNAL(clicked()), this, SLOT(rcv_playRequest()) );

    stepBackButton = new QToolButton(this);
    stepBackButton->setIcon( style()->standardIcon(QStyle::SP_MediaSeekBackward) );
    connect( stepBackButton, SIGNAL(clicked()), this, SLOT(rcv_stepBackwardRequest()) );

    stepForwardButton = new QToolButton(this);
    stepForwardButton->setIcon( style()->standardIcon(QStyle::SP_MediaSeekForward) );
    connect( stepForwardButton, SIGNAL(clicked()), this, SLOT(rcv_stepForwardRequest()) );

//    positionSlider = new QSlider( this );

//    framesPerSecondLabel = new QLabel(tr("Playback Speed"));
//    framesPerSecondSlider = new QSlider();

//    zoomOutButton = new QToolButton(this);
//    zoomOutButton->setText("-");
//    zoomInButton = new QToolButton(this);
//    zoomOutButton->setText("+");

    QBoxLayout *playbackControlsLayout = new QVBoxLayout;
    playbackControlsLayout->addWidget( positionSlider );
    QBoxLayout *playbackButtonsLayout = new QHBoxLayout;
    playbackButtonsLayout->addWidget( playButton );
    playbackButtonsLayout->addWidget( stepBackButton );
    playbackButtonsLayout->addWidget( stepForwardButton );
//    playbackButtonsLayout->addWidget( framesPerSecondSlider );
//    playbackButtonsLayout->addWidget( zoomOutButton );
//    playbackButtonsLayout->addWidget( zoomInButton );
    playbackControlsLayout->addLayout( playbackButtonsLayout );

    setLayout(playbackControlsLayout);
}

void PlaybackParameters::rcv_playRequest() {
    // send play or pause request
    emit send_playRequest( !playing );
}

void PlaybackParameters::rcv_startPauseRequest() {
    pausing = true;
    emit send_playRequest( false );
}

void PlaybackParameters::rcv_endPauseRequest() {
    pausing = false;
    emit send_playRequest( pausing ); // resume previous status
}

void PlaybackParameters::rcv_playingInfo( bool playing ) {
    this->playing = playing;
    QIcon playPause[2] = {
        style()->standardIcon(QStyle::SP_MediaPlay),
        style()->standardIcon(QStyle::SP_MediaPause)
    };
    playButton->setIcon( playPause[ size_t(playing) ] );
}

void PlaybackParameters::rcv_lengthInfo( int length ) {
    this->length = length;
    positionSlider->setRange( 0, length );
    positionSlider->setSingleStep(1);
    positionSlider->update();
    std::cout<<"setRange: 0,"<<length<<"\n";
}

void PlaybackParameters::rcv_positionInfo( int position ) {
    positionSlider->setValue( position );
    this->position = position;
}

void PlaybackParameters::rcv_scrubRequest( int value ) {
    std::cout<<"scrub:"<<value<<"\n";
    position = value;
    emit send_scrubRequest( value );
}

void PlaybackParameters::rcv_stepBackRequest() {
    emit send_stepRequest( true );
}

void PlaybackParameters::rcv_stepForwardRequest() {
    emit send_stepRequest( false );
}
