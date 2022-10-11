
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Siamang-vocal-variation

<!-- badges: start -->
<!-- badges: end -->

This is the repository for D’Agostino et al. (in review).

``` r
us2wav <- tuneR::readWave("/Users/denaclink/Desktop/RStudio Projects/Siamang-vocal-variation/Sound Files/US2.wav")
us2wav <- seewave::cutw(us2wav,from=2.5,to=8,output = 'Wave')

phonTools::spectrogram(us2wav@left,windowlength = 35,fs =44100, 
                       maxfreq = 1400)
```

![](README_files/figure-gfm/unnamed-chunk-2-1.png)<!-- -->

``` r
greatcallwav <- tuneR::readWave("/Users/denaclink/Desktop/RStudio Projects/Siamang-vocal-variation/Sound Files/GreatCall.wav")
greatcallwav <- seewave::cutw(greatcallwav,from=2.5,to=8,output = 'Wave')

phonTools::spectrogram(greatcallwav@left,windowlength = 35,fs =44100, 
                       maxfreq = 1400)
```

![](README_files/figure-gfm/unnamed-chunk-2-2.png)<!-- -->
