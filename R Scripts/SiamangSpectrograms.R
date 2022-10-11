us2wav <- tuneR::readWave("/Users/denaclink/Desktop/RStudio Projects/Siamang-analysis/Sound Files/US2.wav")
us2wav <- seewave::cutw(us2wav,from=2.5,to=8,output = 'Wave')

phonTools::spectrogram(us2wav@left,windowlength = 35,fs =44100, 
                       maxfreq = 1400)

greatcallwav <- tuneR::readWave("/Users/denaclink/Desktop/RStudio Projects/Siamang-analysis/Sound Files/GreatCall.wav")
greatcallwav <- seewave::cutw(greatcallwav,from=2.5,to=8,output = 'Wave')

phonTools::spectrogram(greatcallwav@left,windowlength = 35,fs =44100, 
                       maxfreq = 1400)
