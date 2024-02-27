
# documentation for **neos.py**

**neos** is a python program used at the neutron reflectometer *Amor* at *SINQ* 
to turn *raw data*
into reduced an *orso* comptible *reflectivity file*. 

raw: *nexus hdf5* format:

- event arrays for detector and monitor events
- arrays for device propertise
- entries for instrument configuration

reduced: *orso reflectivity* format:

- header with information about

  - data origin
  - measurement conditions
  - reduction steps

- array with basic or expanded reflectivity data (in the simplest version: the *reflectivity curve*)

## environment

**neos** (version 2.0 and later) was developen with python3.9. 

The following (non-trivial) python modules are required:

- `numpy`
- `h5py`
- `orsopy`

## usage

**neos** is using command line arguments to

- find the raw data
- overwrite default values
- define the parameter range for reduction
- define reduction steps
- define the output path and name

More detailed explanations for some arguments are given in the next section.

### output for `python eos.py -h`:

``` man
usage: neos.py [-h] [-d DATAPATH] [-Y YEAR]
               [-n FILEIDENTIFIER [FILEIDENTIFIER ...]]
               [-r NORMALISATIONFILEIDENTIFIER [NORMALISATIONFILEIDENTIFIER ...]]
               [-sub SUBTRACT] [-o OUTPUTNAME]
               [-of OUTPUTFORMAT [OUTPUTFORMAT ...]]
               [--offSpecular OFFSPECULAR] [-a QRESOLUTION]
               [-ts TIMESLIZE [TIMESLIZE ...]] [-s SCALE [SCALE ...]]
               [-S AUTOSCALE AUTOSCALE] [-l LAMBDARANGE LAMBDARANGE]
               [-t THETARANGE THETARANGE] [-T THETARANGER THETARANGER]
               [-y YRANGE YRANGE] [-q QZRANGE QZRANGE] [-cs CHOPPERSPEED]
               [-cp CHOPPERPHASE] [-co CHOPPERPHASEOFFSET] [-m MUOFFSET]
               [-mu MU] [-nu NU] [-sm SAMPLEMODEL]

eos reads data from (one or several) raw file(s) of the .hdf format, performs
various corrections, conversations and projections and exports the resulting
reflectivity in an orso-compatible format.

optional arguments:
  -h, --help            show this help message and exit

input data:
  -d DATAPATH, --dataPath DATAPATH
                        relative path to directory with .hdf files
  -Y YEAR, --year YEAR  year the measurement was performed
  -n FILEIDENTIFIER [FILEIDENTIFIER ...], --fileIdentifier FILEIDENTIFIER [FILEIDENTIFIER ...]
                        file number(s) or offset (if negative)
  -r NORMALISATIONFILEIDENTIFIER [NORMALISATIONFILEIDENTIFIER ...], --normalisationFileIdentifier NORMALISATIONFILEIDENTIFIER [NORMALISATIONFILEIDENTIFIER ...]
                        file number(s) of normalisation measurement
  -sub SUBTRACT, --subtract SUBTRACT
                        R(q_z) curve to be subtracted (in .Rqz.ort format

output:
  -o OUTPUTNAME, --outputName OUTPUTNAME
                        output file name (withot suffix)
  -of OUTPUTFORMAT [OUTPUTFORMAT ...], --outputFormat OUTPUTFORMAT [OUTPUTFORMAT ...]
  --offSpecular OFFSPECULAR
  -a QRESOLUTION, --qResolution QRESOLUTION
                        q_z resolution
  -ts TIMESLIZE [TIMESLIZE ...], --timeSlize TIMESLIZE [TIMESLIZE ...]
                        time slizing <interval> ,[<start> [,stop]]
  -s SCALE [SCALE ...], --scale SCALE [SCALE ...]
                        scaling factor for R(q_z)
  -S AUTOSCALE AUTOSCALE, --autoscale AUTOSCALE AUTOSCALE
                        scale to 1 in the given q_z range

masks:
  -l LAMBDARANGE LAMBDARANGE, --lambdaRange LAMBDARANGE LAMBDARANGE
                        wavelength range
  -t THETARANGE THETARANGE, --thetaRange THETARANGE THETARANGE
                        absolute theta range
  -T THETARANGER THETARANGER, --thetaRangeR THETARANGER THETARANGER
                        relative theta range
  -y YRANGE YRANGE, --yRange YRANGE YRANGE
                        detector y range
  -q QZRANGE QZRANGE, --qzRange QZRANGE QZRANGE
                        q_z range

overwrite:
  -cs CHOPPERSPEED, --chopperSpeed CHOPPERSPEED
                        chopper speed in rpm
  -cp CHOPPERPHASE, --chopperPhase CHOPPERPHASE
                        chopper phase
  -co CHOPPERPHASEOFFSET, --chopperPhaseOffset CHOPPERPHASEOFFSET
                        phase offset between chopper opening and trigger pulse
  -m MUOFFSET, --muOffset MUOFFSET
                        mu offset
  -mu MU, --mu MU       value of mu
  -nu NU, --nu NU       value of nu
  -sm SAMPLEMODEL, --sampleModel SAMPLEMODEL
                        1-line orso sample model description
```

## read data file(s)

### absolute minimum

purposes:

- fast access to human-readable meta data in the output header
- get an idea about $q_z$ range and statistics

actions:

- read in one raw data file 
- convert the event stream into an $I(\lambda, \alpha_f)$ map
- project this map onto $q_z$ to give an $I(q_z)$ curve
- write this curve in *orso* format to disk

> example:
> 
> `> python neos.py -n 456 -o foo`  
> looks for the file `amor<year>n000456.hdf` in one of the default locations
> (`./`, `./raw/`, `../raw`, local raw data directory on Amor) and
> writes the output to `foo.Rqz.ort`.

### normalisation

purposes:

- fast access to human-readable meta data in the output header
- get a reduced and (partially) corrected reflectivity curve

actions:

- read in raw data file(s) and raw normalisation file(s)
- convert the normalisation measurement into a $N(\lambda, z_\mathrm{detector})$
  map containing information about guide and detector efficiencies,
  illuminated detector area and incoming intensity.
- convert the event stream into an $I(\lambda, \alpha_f)$ map
- normalisation: $R(\lambda, \alpha_f)_{la} = I(\lambda, \alpha_f)_{la} / N(\lambda, z_\mathrm{detector})_{la}$.
- project this map onto $q_z$ to give a $R(q_z)$ curve (not necessarily scaled)
- write this curve in *orso* format to disk

> example:
> 
> `> python neos.py -n 456 -r 123 -o foo`  
> looks for the files `amor<year>n000456.hdf` (reflectivity) and `amor<year>n000123.hdf`
> (normalisation) in one of the default locations
> (`./`, `./raw/`, `../raw`, local raw data directory on Amor) and
> writes the output to `foo.Rqz.ort`.

### read multiple files

- **for the same instrument parameter set (same $\mu$)**

  The arguments of the keys `-n` and `-r` have the general form  
  `<start1>[-<end1>[:<increment1]][,<start2>[-<end2>[:<increment2]],...]`  
  Each number range is defined by a star value, an optional stop value and an
  optional increment. Various ranges are separated by a `,`.

  > example:  
  >
  > `20-25:2,28-30,40` resolves into the list `[20, 22, 24, 28, 29, 30, 40]`

  action:  

  Effectively, the event streams found in the the various files are merged and 
  processed together.

- **for different parameter sets, or to prevent merging**

  The key `-n` accepts more than one argument of the type defined above. The
  (set of) data file(s) related to one argument are merged and give one
  reflectivity curve (one `data_set`) in the output file. The reflectivity
  curves for more than one argument are separated in the output file
  by the separator `# data-set: <identifier>`.

  > example:  
  >
  > `> python neos.py -n 20,21 30 -r 123 -o foo`   
  > results in two reflectivity curves, the first made from files #20 and #21, 
  > the second from file #30. Both are saved in `foo.Rqz.ort`. 

  warning:  
  
  `-r` does accept only one argument! 

### misc.

**year**: The raw file name is created using the file number and the actual year. In case 
the data to be processed were recorded in a previous year, this must be 
explicitely stated with   
`-Y <year>`.

**path**: The default location for the output (and for starting the search for the input files)
iis the present working directory. This can be altered by using the argument   
`-d <path>`.

**subtract $q_z$-dependent curve**:
It is possible to provide a $R(q_z)$ curve in `.Rqz.ort` format to be subtracted 
from the reduced data. E.g. to emphasize the high-$q$ region on a linear scale, or
to illustrate changes in a series of measurements. The argument is   
`-sub <filename>`.

## output options

### other formats

Besides the default *orso* format `.Rqz.ort`, there is the option to
write the $R(\lambda, \alpha_f)$ array and the related input, normalisation and
mask arrays. This output can help with the sample alignment or the readjustment
of parameters (see below), or it can be used for debugging the data processing or 
instrument operation. The suffix of this output is `.Rlt.ort`

The format is chosen by using one or several of the arguments

- `Rqz.ort`, `Rlt.ort`, `Rqz.orb`, `Rlt.orb`
- `Rqz` (= `Rqz.ort` and `Rqz.orb`)
- `Rlt` (= `Rlt.ort` and `Rlt.orb`) 
- `ort` (= `Rqz.ort` and `Rlt.ort`)
- `orb` (= `Rqz.orb` and `Rlt.orb`) 

where `.orb` will be the future nexus comptibe output format.

### $q_z$ binning

The $R(\lambda, \alpha_f)$ arrays are projected onto a $q_z$ grid with bin boundaries
defined by:

$q_{z\,i} \in [0,\, c,\, 2c,\, 3c,\, \dots \hat\imath c] \qquad \forall \quad q_z \le q_\mathrm{base}$,\quad $\hat\imath = q_\mathrm{base} / c$

$q_{z\,\hat\imath+j} \in [q_\mathrm{base} \cdot (1+c)^1, q_\mathrm{base} \cdot (1+c)^2, \dots q_\mathrm{base} \cdot (1+c)^j \dots \qquad \forall \quad  q_z  >q_\mathrm{base}$ 

$q_\mathrm{base} = 0.1\,\mbox{\AA}^{-1}$ is fixed for the moment.
The *resolution* $c$ can be chosen with `-a` among the values
[0.005, 0.01, 0.02, 0.025, 0.04, 0.05, 0.1, 1]
(this is restricted to ensure a *smooth* transition between the 
linear and exponential regions).

### intensity scaling

The argument `-s <value>` leads to a multiplication of all $R(q_z)$ curves
with `<value>`. This is useful for one curve, only, or in combination with the
`-S` argument.

$R(q_z)$ of the first reflectivity curve can be scaled to 1 in the $q_z$ interval define 
by `-S <start> <stop>`. The following $R(q_z)$ curves are then scaled to match the
respective previous one in the overlapping $q_z$ range. 

### time-slizing

One (combined) data set can be chopped in slizes with the argument  
`-ts <interval> [<start> [stop>]]`
where `<interval>` is the time interval length in seconds. The chopping starts
`<start>` seconds after the start of the measurement (default: 0) and ends at  
`<stop>` seconds (default: end of the measurement).

All the resulting $R(q_z)$ curves are stored in one file, one after the other. An additional
column is added with the start time of the respective slize. This enables fast plotting.

> example:  
> 
> `python -n 20-22 -r 123 -ts 60 1200 4000 -f foo`  
> The event streams of the measurements #20, #21 and #22 are merged. All events before
> $t = 1200\,\mathrm{s}$ with respect to the start of meausrement #20 are discarded.
> Then until $t = 4020\,\mathrm{s}$ (or the end of measurement #22) a $R(q_z)$ curve is generated
> for each $60\,\mathrm{s}$ interval.  

## masksing

The specularely reflected intensity illuminated the detector only in a
limited region. To reduce noise, the **detector region of interest** is reduced by default
to the inner horizontal channels of the detector and the vertical channels 
corresponding to the divergence of the beam.
These values can be overwirtten by using the keys `-y` and `-t` for absolute 
finite angles or `-T` for the angular distance from the detector center.

The **wavelength band** can be limited using `-l`, where the arguments are the 
waveengths with unit `angstrom`.

The **$\mathbf{q_z}$ range** can de limited using `-q`, where the arguments are the momentum 
transfer with unit `1/angstrom`.

## overwrite parameters from the nexus file

purposes:

- debugging
- correctiion of wrong entries (due to communication problems)
- take into account misalignments

```
  -cs CHOPPERSPEED, --chopperSpeed CHOPPERSPEED
                        chopper speed in rpm
  -cp CHOPPERPHASE, --chopperPhase CHOPPERPHASE
                        chopper phase
  -co CHOPPERPHASEOFFSET, --chopperPhaseOffset CHOPPERPHASEOFFSET
                        phase offset between chopper opening and trigger pulse
  -m MUOFFSET, --muOffset MUOFFSET
                        mu offset
  -mu MU, --mu MU       value of mu
  -nu NU, --nu NU       value of nu
  -sm SAMPLEMODEL, --sampleModel SAMPLEMODEL
                        1-line orso sample model description
```
