# Data

## `astro`
This folder contains basic properties of different astrophysical objects required for nuSQuIDS calculations, specifically the sun and the earth density profiles. 
### `sun`
Solar models by John Bahcall which were obtained from `http://www.sns.ias.edu/~jnb/SNdata/Export/BS2005/`. 
The most relevant columns are the
1. normalized radius, ranging from 0 to 1. 
The rest of the columns contain different information
2. radius in km
4. Density in g/cm3
7. hydrogen fraction `rH`, which can be related to the electron fraction (`ye = (0.5)(1 + rH)`)
Other columns are referring to differnt elements in the sun, but current nuSQuIDS does not utilize this information and only assumes an isoscalar sun.

### `earth`
There are two options for the earth model. 
`EARTH_MODEL_PREM.dat` is the default one with just density and electron fraction as a function of the normalized radius.
`EARTH_MODEL_PREM_wIso.dat` in addition includes the isotopic fractions of several nuclei in the earth. 

## `xsections`