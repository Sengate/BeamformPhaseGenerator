/* Created on 2022-03-23
 * Author: Onkabetse Sengate
 * About: utils tool for phase calculations
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#ifndef M_2PI
#define M_2PI 6.283185307179586476925286766559005
#endif

#define D2R (M_PI*2.0/360.0)
#define R2D (360.0/(M_PI*2.0))



void RaDec_to_AltAz(float ra , float dec,float lst, float telescope_lat, float *az, float *alt){
    
     //Convert equatorial coordinates (right ascension - decliantio)
     //to Harizontal coordinates (azimuth-altitude)
     /**
     ** Inputs:
     **         ra                  right accesion       [degrees]
     **         dec                 declination          [degrees]
     **         lst                 local sidereal time  [seconds]
     **         telescope_lat       Telescope latitude   [degrees]
     **
     ** Return:
     **         *az                 Azimuths             [degrees]
     **         *alt                Altitudes            [degrees]
     **/
    
    float ha = lst - ra;
    if (ha > 12.0)
        ha = ha - 24.0;
    else if (ha < -12.0)
        ha = ha + 24.0;
    ha = ha * 15.0;
    
    float DEC, HA, LAT, x_hor, y_hor, z_hor;
    //convert to radians
    DEC = dec * D2R;
    HA = ha * D2R;
    LAT = telescope_lat * D2R;
    
    //Horizontal coordinates
    x_hor = cos(HA) * cos(DEC) * sin(LAT) -  sin(DEC) * cos(LAT);
    y_hor =sin(HA) * cos(DEC);
    z_hor = cos(HA) * cos(DEC) * cos(LAT) + sin(DEC) * sin (LAT);
    float r_az = atan2f(y_hor, x_hor) * R2D + 180.;
    float r_alt = asinf(z_hor) * R2D;
    
    *az = r_az;
    *alt = r_alt;
    if (*az >= 360) *az -= 360;
    
}
    


double Jd0(int year, int month, int day)
//Date to Julian days
{
    double JD;
    double AA,BB;
    
    if (month < 3.0){
        month += 12.0;
        year -= 1.0;
    }
    AA = (double)((int)(year/100.0));
    BB = 2 - AA + (double)((int)(AA/4.0));
    
    JD = BB + (double) ((int)(365.25 * year)) +
        (double)((int)(30.6001 * (month + 1))) + day + 1720994.5;
    
    return JD;
    
}


double jd2(struct timespec time_now){
    //Date to Jd
    //Take in timespec
    
    struct tm* timeinfo;
    timeinfo = localtime(&time_now.tv_sec);
    unsigned int year = timeinfo->tm_year + 1900;
    unsigned int month = timeinfo->tm_mon + 1;
    unsigned int day = timeinfo->tm_mday;
    
    return Jd0(year, month, day);
    
}

float get_lst(struct timespec time_now, float tel_long){
    
    //Convert local time to Local Sidereal time for a given
    //Observation longitude
    /**
     ** Inputs:
     **         struct time_now
     **         tel_long             Observation Longitude [degrees]
     **
     ** Return:
     **         LST                   Local sidereal time  [degrees]
     **/
    struct tm* timeinfo;
    timeinfo = localtime(&time_now.tv_sec);
    unsigned int year = timeinfo->tm_year + 1900;
    unsigned int month = timeinfo->tm_mon + 1;
    unsigned int day = timeinfo->tm_mday;
    
    double jd = Jd0(year, month, day);
    
    double T = (jd - 2451545.0)/36525.0;
    
    // Works if time after year 2000, otherwise T is -ve and might break
    double T0 = fmod((6.697374558 + (2400.051336 * T) + (0.000025862 * T * T)), 24.);
    
    double UT = (timeinfo->tm_hour) + (timeinfo->tm_min / 60.)
    + (timeinfo->tm_sec + time_now.tv_nsec / 1.e9) / 3600.;
    
    double GST = fmod((T0 + UT * 1.002737909), 24.);
    double LST = GST + tel_long / 15.;
    while (LST < 0) {
        LST = LST + 24;
    }
    LST = fmod(LST, 24);
    
    return LST;
    
}




