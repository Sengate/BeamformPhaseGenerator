/*##########################################################################
 *
 *###########################################################################
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
//#include "Coordinates.h"
#include "phase_generator.h"

#define phase_index(ichan,jbeam,kant,nBeams,nAnts) ((nBeams * ichan + jbeam)*nAnts + kant)

#define SEC2HR 1/3600.0;
#define Deg2HR 1/15.0;
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define Deg2Rad M_PI/180.0
#define Rad2Deg 180.0/M_PI
#define C_SPEED 299792458.0

/*-----------HIRAX PARAMETERS------------*/
//Using HARTRAO location
#define HIRAX_LONGITUDE 27.6853931
#define HIRAX_LATITUDE  -25.8897515
#define num_CHANNELS 10
#define num_BEAMS 3
#define num_ANTENNAS 256




void _8bit_quantiser( float* fphases, int data_size, char *quat_phases){
    int i, j;
    float num_levels = (float)(powf(2.0, 8)-1.0);
    
    float fmax = fphases[0];
    float fmin = fphases[0];
    for(i=1; i<data_size; i++)
       {if(fmin>fphases[i])
        fmin=fphases[i];
        if(fmax<fphases[i])
            fmax=fphases[i];}
    float step_size = (fmax-fmin)/num_levels;
    
    for (j=0; j<data_size; j++){
        quat_phases[i] = roundf(fphases[i]/step_size);
        
    }
}

//=========================================================================================
//                             READ  THE INPUT FILES
//=========================================================================================
Parameters *read_input_files(const char *ants_file, const char *beams_file){
    
    //Memory allocation
    struct Parameters *param;
    param = Parameters_make_zeros(num_CHANNELS, num_BEAMS, num_ANTENNAS);
    
    //====================================================
    //    READ ANTENNA FILE
    //===================================================
    FILE *ants_pos; //= fopen(ants_positions,"r");
    if ((ants_pos = fopen(ants_file, "r")) == NULL)
    {
        printf("Error : Unable to open Array positions file\n");
        exit(EXIT_FAILURE);
    }
    float x,y,z;
    
    int count =0;
    while ( ( fscanf(ants_pos, "%f %f %f", &x, &y, &z) ) == 3)
    {
        param->EW_antennas[count] = x;
        param->NS_antennas[count] = y;
        param->H_antennas[count] = z;
        
        count++;
    }
    //printf("%f", param->EW_antennas);
    if ( count != num_ANTENNAS)
        
    {
        printf("Number of Antennas in file %d does not match expected number of antennas %d \n", count, num_ANTENNAS);
    }
    
    //====================================================
    //    READ BEAMS FILE
    //===================================================
    FILE *beams_dir;
    if ((beams_dir= fopen(beams_file, "r"))== NULL)
    {
        printf("Error : Unable to open Beams positions file\n");
        exit(EXIT_FAILURE);
    }
    float RAs, DECs;
    
    int count1=0;
    
    while ( ( fscanf(beams_dir, "%f %f", &RAs, &DECs) ) == 2)
    {
        param->RAs[count1] = RAs;
        param->DECs[count1] = DECs;
        count1++;
    }
    if ( count1 != num_BEAMS)
    {
        printf(" Number of beams in file %d does not match number of expected beams %d \n", count1, num_BEAMS);
    }
    
    //param->num_freq = num_CHANNELS;
    float MAX_FREQ = 810.0;
    float MIN_FREQ = 400.0;
    float chn_res = (MAX_FREQ - MIN_FREQ)/num_CHANNELS;
    for (int i = 0; i <  num_CHANNELS; i++) {
        param->frequencies[i] = (MIN_FREQ + i*chn_res) * 10e6;}
    
    return param;
    Parameters_destroy(param);
}

//======================================================================================
// Compute the Julian day (JD)
//=======================================================================================
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
//     Calculate LST from the local time n and telescope longitude-LST in hours
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
    
    //printf("Current Time is %02d:%02d:%lf \n", hours, minutes, seconds+ time_nsec);
    /*------------------------------------------------*/
    float JD = Julian_Day(year, month, day);
    
   // printf("Julian day is %lf \n", JD);
    
    double T = (JD - 2451545.0)/36525.0;
    //gmst in hours
    double gmst =  (24110.54841 + (8640184.812866*T) + (0.093104*T*T - 0.0000063*T*T*T))*SEC2HR;
    gmst = fmod(gmst, 24.);
    //get current time in hours
    double UT = hours + (minutes/60.) + (seconds + time_nsec)*SEC2HR;
    double GSMT = fmod((gmst + UT * 1.002737909), 24.);
    //get the local sidereal time
    double lon = longitude*Deg2HR
    //printf("Longitude is %lf \n", longitude);
    double LST = GSMT + longitude*Deg2HR;
    while (LST < 0) {
        LST = LST + 24;
    }
    LST = fmod(LST, 24);
    
    printf("LST =: %f \n", LST);
    
    return LST;
}

//======================================================================================
//          Convert RA and DEC to Azimuth and Elevation
//======================================================================================
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

/*=========================================================================================
 * Calculate float geometric delays between antennas
 *=========================================================================================
 *------------------
 *Parameters in
 *-----------------
 *Struct param: Structure poluted with frequencies, antenna positions and beam pointings
 *Struct location: Telescope location information (longitude, latitude and height)
 *struct time_now:
 *----------------------
 * Ouput
 *----------------------
 */
complex_char *calculate_quantised_phases(const Parameters *param, struct telescope *location, struct timespec time_now){

    int ichan, jbeam, kant;
    double az, el;
    double *azimuths = (double*)malloc(sizeof(double) * param->num_beams);
    double *elevations = (double*)malloc(sizeof(double) * param->num_beams);
    //get the LST
    double lst = LST(time_now, location->longitude);
    
    complex_char *out_phases;
    out_phases = char_make_zeros(param->num_freq,param->num_beams,param->num_ants);

    //generate phase matrix
    //float *phase_mat = float_phases(param->num_freq, param->num_beams, param->num_ants);
    for (ichan=0; ichan < param->num_freq; ichan++){
        float constant = 2*M_PI * param->frequencies[ichan]/C_SPEED;
        for (jbeam=0; jbeam < param->num_beams; jbeam++){
            
            //get the beams azimuths and altitudes given the LST, location latitude and ....beam pointings
            RaDec2Altaz(param->RAs[jbeam], param->DECs[jbeam], lst, location->latitude, &azimuths[jbeam], &elevations[jbeam]);
            az = azimuths[jbeam];
            el = elevations[jbeam];
            
            for (kant =0; kant <  param->num_ants; kant++){
                //EW_antennas, NS_antennas: antennas locations in the East and North directions
                float EW_projection = param->EW_antennas[kant] * cos(el) * sin(az);
                float NS_projection =  param->NS_antennas[kant] * cos(el) * cos(az);
                
                const int out_index = phase_index(ichan,jbeam,kant,param->num_beams,param->num_ants);
                float geometric_delay = constant * ( EW_projection + NS_projection);
                float real_phases = cos(geometric_delay);
                float imag_phases = -sin(geometric_delay);
                
                //Quantise the phases
                _8bit_quantiser(&real_phases, out_index, out_phases->real);
                _8bit_quantiser(&imag_phases, out_index, out_phases->imag);
            }
        }
    }
    
    return out_phases;
    free(azimuths);
    free(elevations);
    
    char_phase_destroy(out_phases);
    //fphase_destroy(phase_mat);
}


/*
*=======================================================================================
*Make Complex phases
*========================================================================================
*/
/*complex_phases *calculate_ComplexPhases(const Parameters *param, struct telescope *location, struct timespec time_now){
    
    int ichan, jbeam, kant;
    double az, alt;
    double *azimuths = (double*)malloc(sizeof(double) * param->num_beams);
    double *altitudes = (double*)malloc(sizeof(double) * param->num_beams);
    //get the LST
    double lst = LST(time_now, location->longitude);
    
    complex_phases *fphases;
    phases = ComplexPhases_make_zeros(num_CHANNELS,num_BEAMS,num_ANTENNAS);
    
    //float *phase_mat = float_phases(param->num_freq, param->num_beams, param->num_ants);
    for (ichan=0; ichan < param->num_freq; ichan++){
        float constant = 2*M_PI * param->frequencies[ichan]/C_SPEED;
        for (jbeam=0; jbeam < param->num_beams; jbeam++){
            
            //get the beams azimuths and altitudes given the LST, location latitude and ....beam pointings
            RaDec2Altaz(param->RAs[jbeam], param->DECs[jbeam], lst, location->latitude, &azimuths[jbeam], &altitudes[jbeam]);
            az = azimuths[jbeam];
            alt = altitudes[jbeam];
            
            for (kant =0; kant <  param->num_ants; kant++){
        
                float EW_projection = param->EW_antennas[kant] * cos(alt) * cos(az);
                float NS_projection =  param->NS_antennas[kant] * sin(alt) * cos(az);
                
                const int out_index = phase_index(ichan,jbeam,kant,param->num_beams,param->num_ants);
                float geometric_delay = constant * ( EW_projection + NS_projection);

                fphases->real[out_index] = cos(geometric_delay[out_index]);
                fphases->imag[out_index] = -sin(geometric_delay[out_index]);
                
    return fphases;
    complexPhase_destroy(fphases);
}*/

int main(){
    float lon = 20.5;
    //Get the LST
    double lst;
    struct timespec tv;
    
    
    timespec_get(&tv, TIME_UTC);
    char buff[100];
    //strftime(buff, sizeof buff, "%D %T", gmtime(&tv.tv_sec));
    //printf("Current time: %s.%09ld UTC\n", buff, tv.tv_nsec);
    //printf("Raw timespec.time_t: %jd\n", (intmax_t)tv.tv_sec);
    //printf("Raw timespec.tv_nsec: %09ld\n", tv.tv_nsec);
    
    lst  = LST(tv, lon);
    
    printf("%f \n", lst);
    
    //telescope
    //struct telescope HIRAX = {HIRAX_LONGITUDE, HIRAX_LATITUDE};
    //struct telescope *TEL;
    //TEL = &HIRAX;
    
    
    //Polute parameters structure
    //Frequencies, antenna Positions and Beams cordinates
    //struct Parameters *pp;
    //pp = read_input_files("HIRAX_Antenna_Positions.txt", "Beams.txt");
    
    //for (int i=0; i<num_BEAMS; i++){
      // printf("%f \n", pp->frequencies[i]);}
    
    //Calculate float phases
    //float *fPHASES = float_phases(num_CHANNELS, num_BEAMS,num_ANTENNAS);
    
        //fPHASES = calculate_GeometricDelays(pp,TEL,tv);
    
    //Calculate complex phases
   // struct complex_phases *comPhases;
    //comPhases = calculate_ComplexPhases(fPHASES);
    
    /*int ichan, jbeam, kant;
    for (ichan =0; ichan<num_CHANNELS; ichan++){
        for (jbeam =0; jbeam<num_BEAMS; jbeam++){
            for (kant =0; kant<num_ANTENNAS; kant++){
                const int out_idx = phase_index(ichan,jbeam,kant,num_BEAMS, num_ANTENNAS);
        
               // printf("%f\n", fPHASES[out_idx]);
                
                printf("(%+1.5f + %1.5fj) ",comPhases->real[out_idx],comPhases->imag[out_idx]);
                
            }}}*/
    
    return 0;
}
    
    
    
    

