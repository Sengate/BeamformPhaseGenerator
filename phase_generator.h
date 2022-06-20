//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>


struct telescope{
    float longitude;
    float latitude;
};

/*-----------------------------------------------------------------------------------------*/
typedef struct Parameters{
    
    //get frequencies
    float *frequencies;
    unsigned int num_freq;
    
    //get antennas info
    float *EW_antennas;
    float *NS_antennas;
    float *H_antennas;
    unsigned int num_ants;
    
    //get beams info
    float *RAs;
    float *DECs;
    unsigned int num_beams;
    
}Parameters;


Parameters *Parameters_make_zeros(unsigned int num_Channels, unsigned int num_Beams, unsigned int num_Ants){
    
    Parameters *param;
    param = (Parameters *)malloc(sizeof(Parameters));
    
    param->num_freq = num_Channels;
    param->num_ants = num_Ants;
    param->num_beams = num_Beams;
    
    //Allocate frequencies
    param->frequencies = (float*) malloc (sizeof(float) * num_Channels);
    memset(param->frequencies, 0x00, sizeof(float) * num_Channels);
    
    //Allocate Antennas Positions
    param->EW_antennas = (float*) malloc (sizeof(float) * num_Ants);
    memset(param->EW_antennas, 0x00, sizeof(float) * num_Ants);
    
    param->NS_antennas = (float*) malloc (sizeof(float) * num_Ants);
    memset(param->NS_antennas, 0x00, sizeof(float) * num_Ants);
    
    param->H_antennas = (float*) malloc (sizeof(float) * num_Ants);
    memset(param->H_antennas, 0x00, sizeof(float) * num_Ants);
    
    //Allocate Beams coordinates
    param->RAs = (float*) malloc (sizeof(float) * num_Beams);
    memset(param->RAs, 0x00, sizeof(float) * num_Beams);
    
    param->DECs = (float*) malloc (sizeof(float) * num_Beams);
    memset(param->DECs, 0x00, sizeof(float) * num_Beams);
    
    return param;
    
    free(param);
    
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
/*-----------------------------------------------------------------------------------------*/



/*-----------------------------------------------------------------------------------------*/

//complex phases struct
typedef struct complex_phases{
    
//unsigned int nChan;
    //unsigned int nBeam;
    //unsigned int nAnt;
    
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
    
    //phases->nChan = nChan;
    //phases->nBeam = nBeam;
    //phases->nAnt = nAnt;
    
    unsigned int dimesion = nChan * nBeam * nAnt;
    
    phases->real = (float *)malloc(sizeof(float) * dimesion);
    phases->imag = (float *)malloc(sizeof(float) * dimesion);
    
    memset(phases->real, 0x00, sizeof(float) * dimesion);
    memset(phases->imag, 0x00, sizeof(float) * dimesion);
    
    return phases;
    
}

void complexPhase_destroy(complex_phases *phases) {
    free((void *) phases->real);
    free((void *) phases->imag);
    free((void *) phases);
}


/*-----------------------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------------------*/
float *float_phases(const unsigned int NCHAN, const unsigned int NBEAMS,
                    const unsigned int NANTS){
    unsigned int dimesion = NCHAN * NBEAMS * NANTS;
    
    float *phase_mat = (float *)malloc(sizeof(float) * dimesion);
    memset(phase_mat, 0x00, sizeof(float) * dimesion);
    return phase_mat;
}

void fphase_destroy(float *float_phases) {
    free((void *) float_phases);
}

/*-----------------------------------------------------------------------------------------*/



/*==========================================================================================*/
Parameters *read_input_files(const char *ants_file, const char *beams_file);
double Julian_Day(int year, int month, int day);
double LST(struct timespec current_time, double longitude);
void RaDec2Altaz( double RA, double Dec, double LST, double lat, double *Az, double *Alt);
float *calculate_GeometricDelays(const Parameters *param, struct telescope *location, struct timespec time_now);
complex_phases *calculate_ComplexPhases(float *float_delays);


