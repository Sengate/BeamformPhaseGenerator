# BeamformPhaseGenerator


## Phase generation Code


struct  : `Parameters *param` read files containing the following;
Beam pointings in degrees (floats): `pointings[num_beams][Ra, Dec]`,
Dish locations in meters (floats): `dish_locations[num_ants][X,Y,Z]`,
List of frequency centers in MHz (float): `frequency_map[num_frequencies]`

Struct `timespec`:
The time corresponding to the current set of data (used to compute the LST) as a timespec: `time`

Struct `telescope`:
The base corrdinate of the instrument in latitude and longitude as floats: `instrument_corrdinates[lat,long]`  


Outpust:
  
  
From phase_generator.c

Function <<lst>>:
Compute the LST given time now an the telescope longitude
Function <<RaDec2Altaz>>:
Convert the beams equatorial coordinates to horizontal coordinate

Function <<calculate_GeometricDelays>>:
Calculate the geometrical delays for each freq, beam and antenna:
`Output-geometrical_delays[num_frequencies, num_beams, num_dishes]`

Function <<calculate_ComplexPhases>>;
 calculate complex phases 
`Output-Complex_phases[num_frequencies, num_beams, num_dishes]`



## Quantization function (To be included):

Inputs- complex phases, nbits
```
Quantization_steps = (2 * PI)/2^bits
Quantized_phase =  Quantization_steps * round(phase/Quantization_steps)
```
Quantize the real and imag separately
`Output-Complex_phases[num_frequencies, num_beams, num_dishes, 4bit/8bit]`
