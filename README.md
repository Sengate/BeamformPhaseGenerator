## Phase generation function

Function inputs:
Beam pointings in degrees (floats): `pointings[num_beams][Ra, Dec]`
Dish locations in meters (floats): `dish_locations[num_dishes][X,Y, Z]`
List of frequency centers in MHz (float): `frequency_map[num_frequencies]`
Gains (complex floats): `gains[num_dishes, num_frequencies, num_dishes, 2 polarizations]` ??????
The time corresponding to the current set of data (used to compute the LST) as a timespec: `time`
The base corrdinate of the instrument in latitude and longitude as floats: `instrument_corrdinates[lat,long]`  All feed positions are relative to this corrdinate.

Function output:  
Phases to apply to beams: `phases[num_frequencies, num_beams, num_dishes, polarization, complex 4-bit ints]`

For HIRAX-256:
num_dishes = 256
num_frequencies = 64/GPU or 8/stream (total 1024)
num_beams = ? (likely around 400)

For CHORD:
num_dishes = 512
num_frequencies = 16/GPU or 2/stream (total 2048)
num_beams = `[20, 96]`

I've included some notes below from Erik on how the beamforming kernel works.  Note he calls the `phases` array the `A` matrix.  I'm not sure why in his notes the `A` matrix doesn't include polarization, I need to check with him on that.  In Erik's current code he transposes the `phases` matrix into a new format before sending to the GPU.  This might be needed by your code as well, but since we haven't decided on the exact format yet, I think you should just generate the `phases` matrix as given above for now. 

Also note, we are assuming 4-bit complex for now.  However I think we will likely want to make that 8-bit instead (or at least have the option to switch between 4 and 8-bit).

## Beamforming kernel (Kendrick, Erik)

Kernel calculates:



```
J[b,f,t,p] = G[f,b] sum[d] A[f,b,d] E[t,f,d,p]
```

### Current layout:

```
E[32768 times, 16 frequencies, 512 dishes, 2 polarizations, 2 complex, 4 int bits]
A[16 frequencies, 96 beams, 512 dishes, 4 int bits]
G[[16 frequencies, 96 beams, 32 float bits]
J[96 beams, 16 frequencies, 32768 times, 2 polarizations, 2 complex, 4 int bits]
```

A is preprocessed on CPU to actual GPU layout where dishes are renumbered:

```
dish_new = [dish bit 8, dish bits 012345, dish bits 67]
```

This is required because the GPU kernel will re-map the E field in the
same way to transpose the polarizations and complex components out of
the way. With a more convenient E field layout, the A field would not
need to have its dishes renumbered.

### Ideal layout:

```
E[32768 times, 16 frequencies, 2 polarizations, 2 complex, 512 dishes, 4 int bits]
A[16 frequencies, 96 beams, 512 dishes, 4 int bits]
G[[16 frequencies, 96 beams, 32 float bits]
J[96 beams, 16 frequencies, 32768 times, 2 polarizations, 2 complex, 4 int bits]
```

The main point is that the 512 dishes are contiguous in memory. A, G,
and J would be unchanged. Technically, it suffices to have 128 dishes
(1 cache line of dishes) contiguous.
