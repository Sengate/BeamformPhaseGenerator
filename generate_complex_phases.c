#include <stdio.h>
#include "complexPhases_make_zeros.h"
#include "Parameters_make_zeros.h"
#include "Coordinates.h"
#include "read_files.h"

//==================================================================================================================
//                              CALCULATE PHASES
//==================================================================================================================

complex_phases *calculate_complex_phases(const Parameters *param, const telescope *location, struct timespec time_now)
{   //@params param:  Contains frequencies, Antenna Positions and TAB Beams positions
    //@params tel:    Contains Telescope information (LONGITUDE and LATITUDE)
    //@params time_now: Contains Local time
    
    //PHASE_OUT FORMAT: COMPXPHASE[nfreq][nbeams][nants]
    complex_phases *ComplexPhases;
    
    int ichan, ibeam, iant;
    float az, alt;
    
    //get the LST
    float LST = LST(time_now, location->longitude);
    
    //memory allocation to phases/memset to 0X00
    complx_phases = ComplexPhases_make_zeros(param->nfreq, param->nbeams, param->nants);
    
    for (ichan=0; ichan < param->nfreq; ichan++){
        for (ibeam=0; ibeam < param->nbeams; ibeam++){
            for (iant =0; iant < param->nants; iant++){
                float omega = M_2PI * param->frequencies[ichan];
                
                //convert ra/dec to az/alt
                RaDec_to_AltAz(param->RAs[ibeam], param->DECs[ibeam], location->latitude, LST, &az, &alt);
                
                float angx = omega * param->EW_antennas[iant] * sin(alt) * cos(az)/C_SPEED;
                float angy = omega * param->NS_antennas[iant] * sin(alt) * cos(az)/C_SPEED;
                
                complx_phases->real[ichan * param->nbeams * param->nants
                                 + ibeam * param->nants + iant] = cos(angx+angy);
                
                complx_phases->imag[ichan * param->nbeams * param->nants
                                 + ibeam * param->nants + iant] = -sin(angx + angy);
                
            }
        }
    }
    complexPhase_destroy(com_phases);
    return com_phases;
    
}


//==========================================================================================================//
//                        QUANTIZE PHASES
//==========================================================================================================//

complex_phases *PhaseQuantization(const complex_phases *phase, int number_bits){
    
    /*Function takes in Complex Phases and Number of Bits
     *Return Quantized Complex Phases
     *Quant_steps = (2*PI/2**b)
     *Q(Phase) = Quant_steps * roundf(Phase/Quant_steps)
     **/
    complex_phases *Quantized_phases;
    //Quantization steps
    float phaseQuantResolution = 2 * M_PI/ (pow(2, number_bits));
    int ichan, iant, ibeam;
    //make zeros
    Quantized_phases = ComplexPhases_make_zeros(phase->nChan, phase->nBeam, phase->nAnt);
    
    for (ichan=0; ichan < phase->nChan; ichan++){
        for (ibeam=0; ibeam < phase->nBeam; ibeam++){
            for (iant =0; iant < phase->nAnt; iant++){
                float Phase_real = phase->real[ichan*phase->nBeam*phase->nAnt + ibeam * phase->nAnt + iant];
                float Phase_imag = hase->imag[ichan*phase->nBeam*phase->nAnt + ibeam * phase->nAnt + iant]
                
                /*------Quantize Phases---------------*/
                Quantized_phases->real[ichan * phase->nBeam * phase->nAnt
                                       + ibeam * phase->nAnt + iant] = phaseQuantResolution * roundf(Phase_real/phaseQuantResolution);
            
                Quantized_phases->imag[ichan * phase->nBeam * phase->nAnt
                                       + ibeam * phase->nAnt + iant] = phaseQuantResolution * roundf(Phase_imag/phaseQuantResolution);
                
            }
        }
    }
    
    complexPhase_destroy(Quantized_phases);
    return Quantized_phases;
    
}

