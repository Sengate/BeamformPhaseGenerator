#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// A structure holding relevant observation information
//can be used to read info from metafits file

typedef struct Parameters{
    
    //get frequencies
    float *frequencies;
    unsigned int nfreq;
    
    //get antennas info
    float *EW_antennas;
    float *NS_antennas;
    float *H_antennas;
    unsigned int nants;
    
    //get beams info
    float *RAs;
    float *DECs;
    unsigned int nbeams;
    
}Parameters;


Parameters *Parameters_make_zeros(unsigned int nBeams, unsigned int nAnts, unsigned int nChan){
    
    Parameters *param;
    param = (Parameters *)malloc(sizeof(Parameters));
    
    param->nfreq = nChan;
    param->nants = nAnts;
    param->nbeams = nBeams;
    
    
    
    //Allocate frequencies
    param->frequencies = (float*) malloc (sizeof(float) * nChan);
    memset(param->frequencies, 0x00, sizeof(float) * nChan);
    
    
    //Allocate Antennas Positions
    param->EW_antennas = (float*) malloc (sizeof(float) * nAnts);
    memset(param->EW_antennas, 0x00, sizeof(float) * nAnts);
    
    param->NS_antennas = (float*) malloc (sizeof(float) * nAnts);
    memset(param->NS_antennas, 0x00, sizeof(float) * nAnts);
    
    param->H_antennas = (float*) malloc (sizeof(float) * nAnts);
    memset(param->H_antennas, 0x00, sizeof(float) * nAnts);
    
    //Allocate Beams coordinates
    param->RAs = (float*) malloc (sizeof(float) * nBeams);
    memset(param->RAs, 0x00, sizeof(float) * nBeams);
    
    param->DECs = (float*) malloc (sizeof(float) * nBeams);
    memset(param->DECs, 0x00, sizeof(float) * nBeams);
    
    return param;
    
}



//Destroy Inputs memory allocation
void Parameters_destroy( Parameters *param){
    free((void *) param->frequencies);
    free((void *) param->RAs);
    free((void *) param->DECs);
    free((void *) param->EW_antennas);
    free((void *) param->NS_antennas);
    free((void *) param->H_antennas);
    
    
}






