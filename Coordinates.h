/* Created on 2022-03-23
 * Author: Onkabetse Sengate
 * About: utils tool for phase calculations
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

//PI
#define SEC2HR 1/3600.0;
#define Deg2HR 1/15.0;

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define Deg2Rad M_PI/180.0
#define Rad2Deg 180./M_PI


//=========================================================================================
//   
//=========================================================================================
double Julian_Day(int year, int month, int day)
/**
 * Func: Compute the Julian day (JD), given Year, Month and Day
 *
 **Param year:      Year
 **Param month:     Month
 **Param day:       Day
 **
 **Return Julian Day: JD
**/
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

//=========================================================================================
//     Calculate LST from local time and telescope longitude-LST in hours
//=========================================================================================
double LST(struct timespec current_time, double longitude){
  
    /*----------- get time information-------------------*/
    clock_gettime(CLOCK_REALTIME, &current_time);
    struct tm *timeinfo = localtime(&current_time.tv_sec);
    unsigned int hours, minutes, seconds, day, month, year;
    year = timeinfo->tm_year + 1900;
    month = timeinfo->tm_mon + 1;
    day = timeinfo->tm_mday;
    hours =  timeinfo->tm_hour;
    minutes =  timeinfo->tm_min;
    seconds = timeinfo->tm_sec;
    double time_nsec = current_time.tv_nsec/1e9;
    
    printf("Time is %02d:%02d:%lf \n", hours, minutes, seconds+ time_nsec);
    /*------------------------------------------------*/
    
    double JD = Julian_Day(year, month, day);
    
    printf("Julian day is %lf \n", JD);
    double T = (JD - 2451545.0)/36525.0;
    //gmst in hours
    double gmst =  (24110.54841 + (8640184.812866*T) + (0.093104*T*T - 0.0000063*T*T*T))/SEC2HR;
    gmst = fmod(gmst, 24.);
    //get current time in hours
    double UT = hours + (minutes/60.) + (seconds + time_nsec)/SEC2HR;
    double GSMT = fmod((gmst + UT * 1.002737909), 24.);
    //get the local sidereal time
    double lon = longitude/Deg2HR
    printf("Longitude is %lf \n", lon);
    
    double LST = GSMT + longitude/Deg2HR;
    while (LST < 0) {
        LST = LST + 24;
    }
    LST = fmod(LST, 24);
    
    return LST;
    
}

//=========================================================================================
//          Convert RA and DEC to Azimuth and Elevation
//=========================================================================================

void RaDec2Altaz( double RA, double Dec, double LST, double lat, double *Az, double *Alt){
    
    double rha, rdec, ralt, raz, rlat;
    
    double HA = LST - RA;
    if (HA > 12.0)
        HA -= 24.0;
    else if (HA < -12.0)
        HA *=15.0;
    
    double sineLat, cosLat, sineDec, cosDec;
    double sineHA, cosHA, sineAlt, cosAlt;
    //convert angels to radians
    rha = HA * Deg2Rad;
    rlat = lat * Deg2Rad;
    rdec = Dec * Deg2Rad;
    
    sineLat = sin(rlat); cosLat = cos(rlat);
    sineDec = sin(rdec); cosDec =cos(rdec);
    sineHA = sin(rha);   cosHA = cos(rha);
    
    double xhor = cosHA * cosDec * sineLat - sineDec * cosLat;
    double yhor = sineHA * cosDec;
    double zhor = cosHA * cosDec * cosLat + sineDec * sineLat;
    
    //get az/alt in degrees
    raz = atan2f(yhor, xhor);
    *Az =  raz * Rad2Deg + 180.0;
    if (*Az >= 360) *Az -= 360;
    ralt = asinf(zhor);
    *Alt = ralt * Rad2Deg;
    
    
}



