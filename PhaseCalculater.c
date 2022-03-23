#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#include <Phases.h>


#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif
#ifndef M_2PI
#define M_2PI 6.283185307179586476925286766559005
#endif

#define D2R (M_PI*2.0/360.0)
#define R2D (360.0/(M_PI*2.0))

#define speed_c 30000000.0



typedef struct Frequencies{
    
    float * frequencies;
    unsigned int nfreq;
}Frequencies;


typedef struct Coordinates{
    
    float * RAs;
    float * Decs;
    float * Altitudes;
    float * Azimuths;
    
    unsigned int nBeams;

}Coordinates;



struct Tel_gps{
    float longitude;
    float latitude;
};

/*---------------------------------------------------------------------------------------------*/

typedef struct Telescope{
    
    float * x_antsPos;
    float * y_antsPos;
    struct Tel_gps telescope_pos;
}Telescope;



/*------------------------------------------------------------------------------------------------*/

/*------------------------------------------------------------------------------------------------*/

Coordinates * coord_zeros( unsigned int nBeams){
    
    Coordinates *cord;
    
    cord = (Coordinates *) malloc(sizeof(Coordinates));
    
    cord->nBeams = nBeams;
    
    cord->RAs = (float *) malloc (sizeof(float) * nBeams);
    memset(cord->RAs, 0x00, sizeof(float) * nBeams);
    
    cord->Decs = (float *) malloc (sizeof(float) * nBeams);
    memset(cord->Decs, 0x00, sizeof(float) * nBeams);
    
    cord->Altitudes = (float *) malloc (sizeof(float) * nBeams);
    memset(cord->Altitudes, 0x00, sizeof(float) * nBeams);
    
    cord->Azimuths = (float *) malloc (sizeof(float) * nBeams);
    memset(cord->Azimuths, 0x00, sizeof(float) * nBeams);
                                  
    return cord;
}

                                  
void Coord_destroy( Coordinates *cord){
    free((void *) cord->RAs);
    free((void *) cord->Decs);
    free((void *) cord->Altitudes);
    free((void *) cord->Azimuths);
    
                                      
    }


/*------------------------------------------------------------------------------------------------*/



typedef struct Phases{
    
    unsigned int nChan;
    unsigned int nBeam;
    unsigned int nAnt;
    
    float *** outPhases;
    
    float *** EWPhases;
    float *** NSPhases;
}Phases;


/*-----------------------------------------------------------------------------------------------*/
/*                            Construct Phases Matrices                                          */
/*-----------------------------------------------------------------------------------------------*/

Phases * phases_contruct_zeros(const unsigned int nChan, const unsigned int nBeam,
                               const unsigned int nAnt)
{
    Phases *phs;
    
    unsigned int iChan;
    unsigned int iBeam;
    unsigned int iAnt;
    
    phs = (Phases *) malloc(sizeof(Phases));
    phs->nChan = nChan;
    phs->nBeam = nBeam;
    phs->nAnt = nAnt;
    
    phs->outPhases = (float***) malloc(sizeof(float**) * phs->nChan);
    phs->EWPhases = (float***) malloc(sizeof(float**) * phs->nChan);
    phs->NSPhases = (float***) malloc(sizeof(float**) * phs->nChan);
    
    
    for (iChan=0; iChan < phs->nChan; iChan++) {
        phs->outPhases[iChan] = (float**)malloc(sizeof(float*) * phs->nBeam);
        for (iBeam=0; iBeam < phs->nBeam; iBeam++)
            
        phs->outPhases[iChan][iBeam] = (float*)malloc(sizeof(float)*phs->nAnt);
        phs->EWPhases[iChan][iBeam] = (float*)malloc(sizeof(float)*phs->nAnt);
        phs->NSPhases[iChan][iBeam] = (float*)malloc(sizeof(float)*phs->nAnt);
        
        memset(phs->outPhases, 0x00, phs->outPhases[0][0][0] * phs->nChan* phs->nBeam * phs->nAnt);
        memset(phs->EWPhases, 0x00, phs->outPhases[0][0][0] * phs->nChan* phs->nBeam * phs->nAnt);
        memset(phs->NSPhases, 0x00, phs->outPhases[0][0][0] * phs->nChan* phs->nBeam * phs->nAnt);
        
    }
    
    return phs;
}


/*-----------------------------------------------------------------------------------------------*/
/*                            Destroy Phases                                                     */
/*-----------------------------------------------------------------------------------------------*/
void Phases_destroy(Phases *phs) {
    
    unsigned int iChan;
    unsigned int iBeam;
    
    for (iChan = 0; iChan < phs->nChan; iChan++)
    {
        for (iBeam = 0; iBeam< phs->nBeam; iBeam++){
            free((void *) phs->outPhases[iChan][iBeam]);
            free((void *) phs->EWPhases[iChan][iBeam]);
            free((void *) phs->NSPhases[iChan][iBeam]);
        }
        free((void *) phs->outPhases[iChan]);
        free((void *) phs->EWPhases[iChan]);
        free((void *) phs->NSPhases[iChan]);
    }
    free((void *) phs->outPhases);
    free((void *) phs->EWPhases);
    free((void *) phs->NSPhases);
    free((void *) phs);
    
}




/*-----------------------------------------------------------------------------------------------*/
/*                            Calculate Phases                                                     */
/*-----------------------------------------------------------------------------------------------*/

Phases *calculatePhases(const Coordinates *coord, const Frequencies *freq, const Telescope *ants, float lst){
    
    Phases *phs;
    
    unsigned int iBeam;
    unsigned int iChan;
    unsigned int iAnt;
    
    phs = phases_contruct_zeros(freq->nfreq ,coord->nBeams, ants->nant);
    
    for (iChan=0; iChan < freq->nfreq; iChan++){
        
        float omega[iChan] = M_2PI * freq->frequencies[iChan];
        
        for (iBeam=0; iBeam < coord->nBeams; iBeam++ ){
            
            //complute_alt_az(coord->RAs[iBeam], coord->Decs[iBeam], lst);
            
            for (iAnt=0; iAnt < ants->nant; iAnt++){
                
                float angx = omega * ant->x_antsPos[iAnt] * sin(coord->alt[iBeam]) * cos(coord->az[iBeam])/speed_c
                
                float angy = omega * ant->y_antsPos[iAnt] * sin(coord->alt[iBeam]) * sin(coord->az[iBeam])/speed_c
                
                //phs->outPhases[iChan][iBeam][iAnt] = (separate the imaginary and reals ???)
                //phs->EWPhases[iChan][iBeam][iAnt]=???
                //phs->NSPhases[iChan][iBeam][iAnt]=????
            }
        }
    }
    
    
    //return phs;
    //Phases_destroy(phs)
    
    
}







