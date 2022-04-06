#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>



#include "Coordinates.h"
//#include "Phases_Inputs.h"
//#include "generate_phases.h"


//HIRAX PARAMETERS
#define NFREQ 1024
#define NBEAMS 400
#define NANTS 256
#define C_SPEED 30000000.0


/*---------------------------------------------------------------------------------------*/
// Strut holding relevant observation information
//used to read file
/*---------------------------------------------------------------------------------------*/

typedef struct Inputs{
    
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
    
} Inputs;


Inputs *inputs_zeros(unsigned int nBeams, unsigned int nAnts, unsigned int nChan){
    
    Inputs *inp;
    inp = (Inputs *)malloc(sizeof(Inputs));
    
    inp->nfreq = nChan;
    inp->nants = nAnts;
    inp->nbeams = nBeams;
    
    
    
    //Allocate frequencies
    inp->frequencies = (float*) malloc (sizeof(float) * nChan);
    memset(inp->frequencies, 0x00, sizeof(float) * nChan);
    
    
    //Allocate Antennas Positions
    inp->EW_antennas = (float*) malloc (sizeof(float) * nAnts);
    memset(inp->EW_antennas, 0x00, sizeof(float) * nAnts);
    
    inp->NS_antennas = (float*) malloc (sizeof(float) * nAnts);
    memset(inp->NS_antennas, 0x00, sizeof(float) * nAnts);
    
    inp->H_antennas = (float*) malloc (sizeof(float) * nAnts);
    memset(inp->H_antennas, 0x00, sizeof(float) * nAnts);
    
    //Allocate Beams coordinates
    inp->RAs = (float*) malloc (sizeof(float) * nBeams);
    memset(inp->RAs, 0x00, sizeof(float) * nBeams);
    
    inp->DECs = (float*) malloc (sizeof(float) * nBeams);
    memset(inp->DECs, 0x00, sizeof(float) * nBeams);
    
    return inp;
    
}



//Destroy Inputs memory allocation
void inputs_destroy( Inputs *inp){
    free((void *) inp->frequencies);
    free((void *) inp->RAs);
    free((void *) inp->DECs);
    free((void *) inp->EW_antennas);
    free((void *) inp->NS_antennas);
    free((void *) inp->H_antennas);
    
    
}
/*---------------------------------------------------------------------------------------*/


/*------------------------------------------------------------------------------------------------------*/
 //                   COMPLEX PHASE STRUCT AND MEMORY ALLOCATION
/*------------------------------------------------------------------------------------------------------*/


typedef struct complex_phases{
    
    unsigned int nChan;
    unsigned int nBeam;
    unsigned int nAnt;
    
    float *real;
    float *imag;
    
}complex_phases;




//make zeros of complex phases
complex_phases * cPhases_make_zeros(const unsigned int nChan, const unsigned int nBeam,
                                    const unsigned int nAnt)
{
    complex_phases *phases;
    
    phases =(complex_phases *) malloc (sizeof(complex_phases));
    
    phases->nChan = nChan;
    phases->nBeam = nBeam;
    phases->nAnt = nAnt;
    
    unsigned int dimesion = phases->nChan * phases->nBeam * phases->nAnt;
    
    phases->real = (float *)malloc(sizeof(float) * dimesion);
    phases->imag = (float *)malloc(sizeof(float) * dimesion);
    
    memset(phases->real, 0x00, sizeof(float) * dimesion);
    memset(phases->imag, 0x00, sizeof(float) * dimesion);
    
    return phases;
    
}



//Destroy complex phases
void complexPhase_destroy(complex_phases *phases) {
    
    free((void *) phases->real);
    free((void *) phases->imag);
    free((void *) phases);
}

/*------------------------------------------------------------------------------------------------------*/
/*                READ FILES FUNCTION
/*------------------------------------------------------------------------------------------------------*/

Inputs *read_files(const char *ants_positions, const char *BeamsCords)
{
    Inputs *inputs;
    
    inputs->nants = NANTS;
    inputs->nbeams = NBEAMS;
    inputs->nfreq = NFREQ;
    
    //Allocate memory
    inputs  = inputs_zeros(NBEAMS, NANTS, NFREQ);
    
    FILE *ants_file; //= fopen(ants_positions,"r");
    
    if ((ants_file = fopen(ants_positions, "r")) == NULL)
    {
        printf("Error : Unable to open Array positions file\n");
        exit(EXIT_FAILURE);
    }
    float x,y,z;
    
    int count =0;
    while ( ( fscanf(ants_file, "%f %f %f", &x, &y, &z) ) == 3)
    {
        inputs->EW_antennas[count] = x;
        inputs->NS_antennas[count] = y;
        inputs->H_antennas[count] = z;
        
        count++;
        if ( count >= inputs->nants)
        {
            break;
        }
        
    }
    
    
    /* -----------------Read beams Positions----------------------*/
    
    inputs->nbeams = NBEAMS;
    
    FILE *beamsCoords;
    int ra, dec;
    
    if ((beamsCoords= fopen(BeamsCords, "r")) == NULL)
    {
        printf("Error : Unable to open Beams coordinates file\n");
        exit(EXIT_FAILURE);
    }
    int count2 =0;
    while ( ( fscanf(beamsCoords, "%d %d", &ra, &dec) ) == 2)
    {
        inputs->RAs[count2] = ra;
        inputs->DECs[count2] = dec;
        count2++;
        if ( count2 >= inputs->nbeams)
        {
            break;
        }
    }
    
    
    /*---------------------Get Frequencies-----------------------------*/
    inputs->nfreq = NFREQ;
    float MAX_FREQ = 810.0;
    float MIN_FREQ = 400.0;
    float chn_res = (MAX_FREQ - MIN_FREQ)/NFREQ;
    for (int i = 0; i < NFREQ; i++) {
        inputs->frequencies[i] = (MIN_FREQ + i*chn_res) * 10e6;}
   
    for (int i = 0; i < NFREQ; i++){
        printf("%f",inputs->frequencies[i] );
    }
    
    //memory delocation
    inputs_destroy(inputs);
    
    return inputs;
    
    
}




/*------------------------------------------------------------------------------------------*/
//                              CALCULATE PHASES
/*------------------------------------------------------------------------------------------*/
complex_phases *calculate_complex_phases(Inputs *inputs, struct timespec time_now, float LATITUDE, float LONGITUDE)
{
    
    //calculate complex phases
    
    complex_phases *com_phases;
    
    
    int ichan, ibeam, iant;
    float az, alt;
    //get the LST
    float LST = get_lst(time_now, LONGITUDE);
    
    //memory allocation to phases/memset to 0X00
    com_phases = cPhases_make_zeros(NFREQ, NBEAMS, NANTS);
    
    for (ichan=0; ichan < com_phases->nChan; ichan++){
        
        
        for (ibeam=0; ibeam < com_phases->nBeam; ibeam++){
            
            
            for (iant =0; iant < com_phases->nAnt; iant++){
                
                float omega = M_2PI * inputs->frequencies[ichan];
                
                //convert ra/dec to az/alt
                RaDec_to_AltAz(inputs->RAs[ibeam], inputs->DECs[ibeam], LATITUDE, LST, &az, &alt);
                
                float angx = omega * inputs->EW_antennas[iant] * sin(alt) * cos(az)/C_SPEED;
                
                float angy = omega * inputs->NS_antennas[iant] * sin(alt) * cos(az)/C_SPEED;
                
                com_phases->real[ichan * com_phases->nBeam * com_phases->nAnt
                                 + ibeam * com_phases->nAnt + iant] = cos(angx+angy);
                
                com_phases->imag[ichan * com_phases->nBeam * com_phases->nAnt
                                 + ibeam * com_phases->nAnt + iant] = -sin(angx + angy);
                
                
            }
        }
    }
    
    complexPhase_destroy(com_phases);
    
    return com_phases;
    
}






int main(){
    
    
    

    
    return 0;
}
