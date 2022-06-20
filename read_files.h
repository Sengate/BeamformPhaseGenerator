
#include "Parameters_make_zeros.h"

#define NUM_BEAMS 10
#define NUM_ANTENNAS 10
#define NUM_FREQ 500



//============================================================================================
//                             READ FILES
//============================================================================================
Parameters *read_input_files(const char *ants_file, const char *beams_file){
    Parameters *param;
    param = Parameters_make_zeros(NUM_BEAMS, NUM_ANTENNAS, NUM_FREQ);
    
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
    if ( count != NUM_ANTENNAS)
        
    {
        printf("Number of Antennas in file %d does not match expected number of antennas %d \n", count, NUM_ANTENNAS);
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
    if ( count1 != NUM_BEAMS)
    {
        printf(" Number of beams in file %d does not match number of expected beams %d \n", count1, NUM_BEAMS);
    }
    
    param->nfreq = NUM_FREQ;
    float MAX_FREQ = 810.0;
    float MIN_FREQ = 400.0;
    float chn_res = (MAX_FREQ - MIN_FREQ)/NUM_FREQ;
    for (int i = 0; i < NUM_FREQ; i++) {
        param->frequencies[i] = (MIN_FREQ + i*chn_res) * 10e6;}
    
    return param;
    Parameters_destroy(param);
    
}


void read_antenna_positions(const char *file_name, const Parameters *param){
    
    FILE *ants_file; //= fopen(ants_positions,"r");
    if ((ants_file = fopen(file_name, "r")) == NULL)
    {
        printf("Error : Unable to open Array positions file\n");
        exit(EXIT_FAILURE);
    }
    float x,y,z;
    int count =0;
    while ( ( fscanf(ants_file, "%f %f %f", &x, &y, &z) ) == 3)
    {
        param->EW_antennas[count] = x;
        param->NS_antennas[count] = y;
        param->H_antennas[count] = z;
        
        count++;
    }
    //printf("%f", param->EW_antennas);
    if ( count != NUM_ANTENNAS)
        
    {
        printf("Number of Antennas in file %d does not match expected number of antennas %d \n", count, NUM_ANTENNAS);
    }
    
}


void read_beam_directions(const char *file_name, const Parameters *param){
    
    FILE *beams_dir;
    if ((beams_dir= fopen(file_name, "r"))== NULL)
    {
        printf("Error : Unable to open Beams positions file\n");
        exit(EXIT_FAILURE);
    }
    float RAs, DECs;
    int count=0;
    while ( ( fscanf(beams_dir, "%f %f", &RAs, &DECs) ) == 2)
    {
        param->RAs[count] = RAs;
        param->DECs[count] = DECs;
        count++;
    }
    if ( count != NUM_BEAMS)
    {
        printf(" Number of beams in file %d does not match number of expected beams %d \n", count, NUM_BEAMS);
    }
    
}


