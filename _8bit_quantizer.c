#include <stdio.h>
#include <math.h>
#include <stdint.h>


void _8bit_Quantizer(float *data, float dataMin, float dataMax, int data_size, uint8_t *_8bit_data){
    
    int i, j;
    int nbit = 8;
    float num_Quant_levels = (float)(powf(2.0, nbit)-1.0);
    
    //find data maximum and min
    for (i=0; i<data_size; i++){
        if (data[i]>dataMax)
            dataMax = data[i];
        if (data[i] < dataMin)
            dataMin = data[i];
    }
    //size of quantization levels
    float step_size = (float)(dataMax - dataMin)/num_Quant_levels;
    
    for (j=0; j<data_size; j++){
        _8bit_data[j] =  (uint8_t) roundf((data[i] - dataMin)/step_size)* step_size + dataMin;
    }
    
}

