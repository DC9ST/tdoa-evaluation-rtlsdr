# Matlab Scripts for Evaluation of a TDOA System based on 3 RTL-SDRs

## Functionality:
Matlab scripts according to my project on transmitter localization with time-difference-of-arrival (TDOA).

<http://www.panoradio-sdr.de/tdoa-transmitter-localization-with-rtl-sdrs/>

[Presentation at the Software Defined Radio Academy, 2017](https://www.youtube.com/watch?v=Km4TU17b05s)

Update 2019:
All algorithm parameters including the positions of the receivers etc. can be specified in the configuration .m file 
The script is flexible and can be adapted to a 3-receiver scenario (with reference TX) anywhere.

### Input:
Three recordings from (the three) receivers done with the librtlsdr-2freq (modified librtlsdr for 2 frequency operation).
First frequency is the reference transmitter, second frequency is the measurement frequency (required length: 3.6e6 (i.e. 1.2e6 samples per frequency), 2 Msps sampling rate).

### Processing:
* Delay analysis of reference signal to synchronize the receptions
* Delay analysis of measurement signal to measure the TDOAs
* Available Options: bandwidth filtering, correlation smoothing, correlation type switching, signal interpolation
* Mulitlateration calculations for localization

### Output:
A html/javascript file with a map showing the receivers, the hyperbolas and the most likely position of the transmitter.



