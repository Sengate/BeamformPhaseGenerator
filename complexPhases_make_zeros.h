#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//complex phases struct
typedef struct complex_phases{
    
    unsigned int nChan;
    unsigned int nBeam;
    unsigned int nAnt;
    
    float *real;
    float *imag;
    
}complex_phases;

//make zeros of complex phases
//Dimensions [nChan, nBeams, nAntenna]
complex_phases * ComplexPhases_make_zeros(const unsigned int nChan, const unsigned int nBeam,
                                    const unsigned int nAnt)
{
    complex_phases *phases;
    
    phases =(complex_phases *) malloc (sizeof(complex_phases));
    
    phases->nChan = nChan;
    phases->nBeam = nBeam;
    phases->nAnt = nAnt;
    
    unsigned int dimesion = nChan * nBeam * nAnt;
    
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
    
//Print phases
void complex_phases_print(const complex_phases *phases){
    
    int i,j,k;
    for (i =0; i<phases->nChan; i++){
        for (j =0; j<phases->nBeam; j++){
            for (k =0; k<phases->nAnt; k++){
            printf("(%+1.5f + %1.5fj) ",phases->real[i*phases->nBeam*phases->nAnt + j * phases->nAnt + k],phases->imag[i*phases->nBeam*phases->nAnt + j * phases->nAnt + k]);
            }
            
        }
        
    }
             printf("\n");
}



