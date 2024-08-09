% User Manual 
  for the
  Neutron Reflectometer Amor
% Jochen Stahn
% 2024-06-06

---

# Users' guide to NICOS

## open NICOS

- *normal* users should be logged on as `amorlnsg` on the instrument control computer `amor.psi.ch`.
  That is the one under the desk in the control cabin.  
  
  > computer: `amor.psi.ch`  
  > username: `amorlnsg`  
  > password: `<YY>lns1` where `<YY>` is the present year

- create a local subdirectory for the actual campagne:

  1. Open a terminal and enter your name. Probably create a new directory with your name.
     You will end up in this directory.
  1. Create a sub-directory: `> mkdir <year>-<month>` (or the like).
  1. `> cd <year>-<month>`

- open the NICOS gui by typing 

  > `> nicos-gui`

  in a terminal window as user `amorlnsg`.
  The gui will start up. Probably you will have to connect to the NICOS program: 
     
  > in the upper right corner type the wheel and select *connect*  
  > user: `user`  
  > password: `<YY>lns1` where `<YY>` is the present year

## identify yourself

Upon starting an measurement campagne, you have to enter the following 
information:

- proposal ID
- title
- user name(s)
- affiliation?

These data are used to create a repository for your data and are
also stored in the meta data section of all your files.

## brief intro for NICOS commands

### check and change device parameters

- everything is a *device* 
- each device has a *name*, in most cases a 2 to 3 letter abbreviation  
  e.g. `som` for the sample tilt (sample omega)
- get the **description** for a device:
  
  > `<device>`

- get the **present value** (position, angle, ...):

  > `<device>()`

- get all the paramters (position, unit, limits, ...):

  > ???

- **move** to a new value (and wait until the command is executed):

  > `maw(<device>, <value> [, <device2>, <value2> ...])`

All these actions are also realised in the gui: klick on the device name in the 
right pannel. A pop-up window opens with all the options available.

### counting and scanning


### batch file creation and execution


## define set-up

...
### devices

`read(<device>)`

: reads the value of a *device* 

`maw(<device>, <value> [, <device2>, <value2> ...])`

: (move and wait) moves a *device* to *value* 


The sample must be positioned with the center of its surface at the 
focal point. At the same time it should be in the center of rotation of
the $\omega$ stage. For this reason there are 2 vertical **translation**s which are
close to parallel: 

`SOZ`
: = *S*ample *O*mega stage *Z* position
: Lift of the $\omega$ stage so that its center is in the FP.

`STZ`
: = *S*ample *T*ranslation *Z* direction
: Lift of the sample on the $\omega$ stage to bring it to the
  center of rotation. 
: *This stage is not available at the moment!*

And finally the sample has to be **tilt**ed by using the $\omega$ and $\chi$ stages.

`MU`
: Tilt of the sample relative to the horizon.
: On amor $\mu$ is used to define and probably to describe the sample
  orientation relative to the lab horizon. Since the beam might be convergent
  on the sample surface it is not the neutron's angle of incidence!
  See also *[coordinate system(s) and nomenclature]*.

`NU`
: Rotation of the detector center around the sample position
relative to the lab horizon. This is a combined movement of
detector lift, tilt, (*x* translation) and also affects 
all other devices behind the sample.

Diaphragms

`D<n><m>`
: = *D*iaphragm number *n*, blade or position *m*
: $n \in \{\mathrm V, 1..4\}$ for $v$irtual source and *n*umber of diaphragm.
: $m \in \{\mathrm{T, B, L, R, H, V, Z}\}$ for *t*op, *b*ottem, *l*eft, *r*ight,
  *h*orizontal, *v*ertical and *z*-position.
: not all options rea available for all devices. 

| name    | description               | *m*       | range or values / mm  |
| :---    | :---                      | :---      | :---:                 |
| `DVV`   | virtual source vertical   |           | 0.5 1 2 3 4 6 8 10    | 
| `DVH`   | virtual source horizontal |           | 2 5 10 12 15 20 25 30 |
| `DMF`   | middle focus              |           | slot 1 .. 5           |
| `D1<m>` | behind Selene guide       | `T B L R` | -40 .. +40            |
| `D2<m>` | before sample             | `T B L R` | -40 .. +40            |
| `D2Z`   | '' lift                   |           | -100 .. +100          |
| `D3<m>` | behind sample             | `T B L R` | -40 .. +40            |
| `D3Z`   | '' lift                   |           | -100 .. +100          |
| `D4<m>` | before detector           | `H V`     | +1 .. +140            |

## perform measurement

The $q_z$ range of one measurement is defined by the wavelength
range  $\lambda \in [3.5, 12.5]$ \AA\ and the spread of the angles 
of incidence $\alpha_i$. 
For specular reflecivity, $\alpha_i$ is deduced from the position (i.e. angle)
of the detector and the position on the detector where a neutron is detected.

`count(<mode>=<preset>)`

: count with *mode* 

  - `t` = time for *preset* seconds
  - `m` = monitor for *preset* monitor counts 

`scan(<device>, <start>, <step>, <np>, <mode>=<preset>)`

: scan the *device* from *start* with *np* steps of width *step* with the
counting time defined by *mode* / *preset*  

`run(<scriptname>)`

: *run* the script *scriptname*, who's path is either relative to 
`Exp.scriptpath` or which has an absolute path. The script language is
python. 


# How to perform a measurement

It is assumed here that NICOS is running, the basic data are provided and the 
general set-up is realised.

## identify sample

Enter a name and a description of the sample. This information will
end up in the meta data section.

    NewSample(<sample name>)

It is also advised to give a short sample description in orso model format:

    Sample.orsomodel='<model>'

e.g.

    Sample.orsomodel='air | 20 ( Ni 5 | Ti 7 ) | SiO2'

## align sample

There are 2 sets of (virtual) devices: the ones for aligning the smaple and the ones
to select a certain $q$ range or instrument set-up.
The **sample alinment devices** are:

`mud`
: for *pitch correction*:
: To compensate the misalignment of the sample surface with respect to the 
: sample table surface along the beam.

`sch`
: for *roll correction*:
: To compensate the misalignment of the sample surface with respect to the 
: sample table surface normal the beam.  

`sz`
: for *height adjustment*:
: To bring the sample surface to the center of rotation = focal point
: of the neutron beam.

It is not necessary - nor is it allowed - to rerdefine any device's zero position!

### height alignment 

In case the sample surface or a mark outside the sample enviromnet indicating
its position are visible, one can use the **surveyor's optical level** mounted 
outside the Amor area at the wall next to the experimental cabin.
This ensures that the sample is in the beam with a vertical offset of a few tenths 
of a mm.

If possible, an **absorber** at a defined distance above (or below) the sample surface
can be used: Scanning the sample though the direct beam allows for detecting this absorber
and moving to the right position.

Using the full divergence and a large virtual source allows for pitch-misalignments
of up to $1^\circ$ and height-offsets of 5 mm to still generate a signal on the detector
for a reflected beam. A good choice for `mu` for searching some signal in this case is
$0.8^\circ$. 

### pitch

If one of the schemes 'simple' or '???' is used, the detector is positioned in a way
that the refelcted beam footprint is centered. This means that the lower edge of
the footprint - which is much easire to detect then the upper edge - is half the divergence
below the detector center. The center is at channel number $z = 224$.

If the full divergenc of $1.5^\circ$ is used, the lower footprint boundary is
at channel number $z = 384$.

### roll

Once *height* and *pitch* are (roughly) aligned, the *roll* can be corrected by looking at the
$I(y, \theta)$ detector image of the reflected beam. The lower boundary of the 
footprint should be exactly horizontal. 

The *roll* is independent of *pitch* and *height*. An iterative alignment process is thus not
necessary.

An accuracy of about $0.1^\circ$ is sufficient.

### fine alignment

The beam cross section at the focal plane is mimimum for the settings

> `div` = $0.05^\circ$  
> `kad` = $0.70^\circ$  
> `dmf` = 1

I.e. *height* and *pitch* alignment with these settings at a moderate $q$ should
give the best results. The iterative procedure is similar to 'conventional'
reflectometers:

> - scan of `sz` to find the *height* maximum
> - adjust `mud` so that the beam footprint if centered at detector
>   channel $z = 224$  

## measurement: ranges, angles and duration



## reduce data

The *raw* data at amor are stored in the nexus hdf format. Besides some meta data
bout the experiment and the set-up, it contains the measurement data in form of an event.
I.e. each neutron has one entry with the time and location (better: voxel-ID) of detection. 

There is the python program code **eos** to reduce the event stream to something like a 
Intensity vs. wavelength and angle map. This process includes a variety of corrections
(instrument geometry, chopper properties, filtering in wavelength and footprint size, gravity, 
etc. If available, a reference measurement (characterising the incident beam, the detector
efficiency and absorption of the sample environment) can be used to create a 
$R(\lambda, \theta)$ map for specular reflectometry, .... 
and finally an $R(q_z)$ curve.

The output format follows the *orso* rules. 

### location of the raw data files:

If not specified differently, **eos* looks for the raw data files in the following directories:

`./`, `./raw/`, `../raw`, `../../raw`,  `/home/amor/data/<year>/<proposalname>`



The present approach is to create a subdirectory for the present experiment

    /home/amorlnsg/<user>/<date_description>

and in there a link to the raw files (if it does not already exist)

    cd /home/amor/data/<year>/<proposalname>
    ln -s /home/amor/data/<year>/<proposalname> raw

### activate the virtual environment for eos

in the terminal in which the data reduction is to be performed type:
 
```bash
source ~/eos/venv/bin/activate.csh
```

### run eos

**eos** is a command line script. This means that all necessary information is
provided as arguments after the program call on the same line.
At least a file (measurement) number and the name of the output file have to be
provided. This then looks like

```bash
python ~/eos/eos.py -f <number> -o <name>
```

All 'missing' information is guessed or collected by eos, where most default values
are extracted form the raw data file (e.g. beam divergence, instrument settings,
sample and detector angles). Almost all relevant values can be overwritten by
an command line option - in most cases this is not recommended!

Useful arguments for data reduction are:

data file(s):

`-f <value>` 
: number(s) of the raw data files  
: the general format is   
: `<number0>-<number1>:<increment0>[,<number2>-<number3>:<increment2>..]`  
: !!!

`-n <value>`
: number of the reference file for normalisation

`-Y`
: year of the data collection (in order to create the correct raw file name.
: The default is the actual year.

filtering / region of interest:

`-l <low> <high>`
: wavelength range in angstroms

`-t <low> <high>`
: theta range

`-y <low> <high>`
: horizontal pixel range on the detector

`-q <low> <high>`
: $q_z$ range

correct angle:

`-m`
: offset to the `mu` stored in the data file

`-mu`
: overwrite the data file entry

scaling:

`-S <low> <high>`
: scale the (weighted) average in the given $q_z$ range to 1

`-s <value>`
: multiply *reflectivity* by *value*; executed after `-S`

formats:

`-of Rqz | Rlt | ort | orb`  
: Defines what is written to the file: 
: `Rqz`: $R(q_z)$, 
: `Rlt`: $R(\lambda, \theta)$ and much more 
: and in which format:
: `ort`: ASCII file
: `orb`: nexus conform binary file

`-h`
: get help on eos and a description of all arguments.



# Instrument description

Amor is a neutron reflectometer with a beam focused to the sample position.

## Parameters and options

- wavelength range
- angular range
- q range
- polarisation efficiency
- angle of incidence on liquid surfaces

- sample environment
  - 1 T electomagnet (with closed cycle refrigerator)
  - 7 T cryomagnet, horizontal field direction
  - furnaces
  - LB trough
  - sheer cell
  - potentiostate / galvanostate


## From physical source to virtual source

The *real neutron source* is the spallation target at SINQ. The fast
neutrons created there by proton capture in lead nuclei are moderated
to room temperature using D~2~O, and further down to cold
neutrons by a liquid H~2~ moderator. This is often referred to as
*cold source*. 

Some of these cold neutrons are guided by a 4.5 m long converging 
neutron guide (with m = 2.5 coating) to the *virtual source* (VS) position, 
just outside the shileding monilith at 6 m from the surface of the
cold moderator.

virtual source

: _
: 2 Boron-Aluminium wheels with horizontal
  and vertical slits of various sizes, respectively.
  The opening (luminous field diaphragm) is centered at the guide
  focal point.
  

Table: Vertical and horizontal slit widths of the virtual source diaphragm.

 | opening | vertical | horizontal  |
 | ---:    | ---:     | ---:        |
 |         |  / mm    |  / mm       |
 |  1      |   0.5    |   2.0       |
 |  2      |   1.0    |   5.0       |
 |  3      |   2.0    |  10.0       |
 |  4      |   3.0    |  12.0       |
 |  5      |   4.0    |  15.0       |
 |  6      |   6.0    |  20.0       |
 |  7      |   8.0    |  25.0       |
 |  8      |  10.0    |  30.0       |

## Selene guide 

The divergent neutron beam emerging from the *virtual source* slits is
collected and focussed to a point (the *mdiddle focus*, MF) 
some 15 m away from the VS by a set of
two elliptically bent mirrors (i.e. a Montel optics), where one reflects horizontally towards
the Aare, and the other vertically upwards. 

The image of the VS at the MF is distorted due to come aberatin. To correct for this
to first order, a second Montel optics follows which reflects towards Berg and downwards.
This results in an image of the VS at the final focal point, some 30 m behind the VS.
This image is distorted spatialy due to imperfections of the Montel optics (surface and alignement)
and it is inhomogeneous in intensity due to an only partially corrected divergence distribution and
due to losses during reflection.

This arrangement of two Montel optics we called **Selene guide**.

The Selene guide is made up of about 500 mm long L-shaped elemets. These are individually
tiltable in vertical direction (*pitch* movement). 6 of these elements are bundeled on a
support beam which rests on two vertically movable sockets. 3 beams form one Montel optics.
All these movements are highly delicat because they might lead to collissions and thus 
demage at places that are not accessable for years - if ever again. Thus no user is
allowd to run these.

The housings of each Montel optics plus the ones of some beam shaping elements and of 
one chopper disk each form a continuous vacuum vessel.

The Montel optics are each some 9 m long and located symmetrically between the adjunct focal
points. This means there are gaps between VS and guide entry and the guide exit to FP of
3 m, and of 6 m between the Montel optics. These gaps are used for conditioning the 
neutron beam.

### characteristic measures

- $a = 7'421\,\mathrm{mm}$
- $b = 129.50\,\mathrm{mm}$
- $c = \sqrt{a^2 - b^2} = 7'420\,\mathrm{mm}$
- $4c = 29'680\,\mathrm{mm}$
- $\Delta \alpha = \arctan \frac{b}{a} \sqrt{\frac{1+\xi}{1-\xi}} 
                 - \arctan \frac{b}{a} \sqrt{\frac{1-\xi}{1+\xi}}$  
  $\approx\, 1.5\, b/a \,\approx\, 1.5^\circ$ for $\xi = 0.6$

## Optics in the bunker

The gap between VS and Montel optics is bridged by an evacuated flight tube.

The large central gap hoists:

polariser / frame overlap filter

: \hfill MF - 2700
: The polariser consists of Fe/Si supermirror coated Si blades which are bent in
  the shape of a logarithmic spiral with the MF as the spiral pole. To enhance 
  the polarisation efficiency, 2 of these sheets are mounted with a distance of a few
  cm. The total arrangement is some 2.1 m long.
: The frame overlap filter (FOF) has a similar design, but only one sheet, coated on the 
  outside (towards the source) with Ni.
: Both are monted one baove the other and switching means a vertical translation of the 
  assembly.

first chopper disk

: (master) \hfill MF - 500
: The chopper disk has 2 openings of 13.5^o^, each, and rotates with a speed of
  1000 rpm creating pulses ever 33.33 s.
: The beam shaping properties are discussed below in context with the second chopper disk.

main shutter

: \hfill MF - 300

beam monitor

: \hfill MF - 250
: The monitor is a fission chamber with sensitivity 10^-?^.

neutron camera

: \hfill MF
: For alignment purposes. By default not in the beam.

middle focus aperture

: \hfill MF
: Wheel with 5 freely configurable slots.
  Presently unused.

second chopper disk

: (slave) \hfill MF + 500
: Geometrically identical to the first chopper disk, but rotating
  in the opposite sense with a phase offset of -13.5^o^.
  This way they form a *blind double chopper* with a resolution
  $\Delta t/t = \Delta \lambda/\lambda = constant$. [@vanWell]
  The center of this set-up, i.e. MF, is the virtual origin in time and space 
  for the time-of-flight encoding of wavelength $\lambda$. 
  

RF spin flipper

: _

and some shielding elements and flight tubes.

These optical elements produce a neutron beam which

 - is convergent to the FP
 - is restricted to a wavelength range of 3 to 9.5 Aa
 - might be polarised with selectable polarisation
 - has a time focus 15 m before the FP (i.e. the pulse they belong to was 
   virtually created 15 m upstream for all divergences due to the 
   equal trajectory lengths in an ellipse). 

## Optical bench

Behind the exit of the second Montel optics the following 
components are / will be / might be installed:

instrument shutter

: _

laser system for sample alignment

: (to come)

diaphragm D1 

: (to come)
: for defining vertical and horizontal divergences
  and for changeing the angle of incidence for a reduced beam  within the
  full divergence of 1.5 deg.
: individual movement of $\pm40$ mm of all 4 blades relative to center

deflector

: (to come, optional)
: to redirect a restricted beam downwards to
  liquid surfaces

diaphragm D2

: (to come)
: to reduce background
: individual movement of $\pm40$ mm of all 4 blades relative to center

sample table

: _

diaphragm D3

: (to come)
: to reduce background
: individual movement of $\pm40$ mm of all 4 blades relative to center

analyser spin flipper

: (future option)

analyser

: (future option)

diaphragm D4

: (to come)
: vertical / horizontal slit with 1 to 160 mm opening, centered

detector

: _
: the active region starts 4'000 mm behind the FP and 18'842 mm
  behind the MF in standard position
: The detector was designed and built by the ESS detector group.
  It is a prototype for the detector for Estia@ESS with a reduced active area
  of $140 \times 160$ mm^2^. The full beam footprint is 
  $110 \times 120$ mm^2^. 

In addition there might be flight tubes and shileding.

## Infrastructure

### Areas

### Media

### EDV


#### Amor network

By default all LAN sockets in the experimental areas (below the granite
beam and against the bunker wall on the upper area) and all
LAN sockets in the cabin below the desk are connected to the Amor
network.

In this network there are also the instrument controll computer
and computers for the detector, the chopper, ....

#### PSI network

By default the 2 LAN sockets in the corner opposite the the doors in the 
cabin are connected to the PSI web.

!!! WLAN: EDUROM, CORP

### PSYS

# Experiments

## coordinate system(s) and nomenclature

The lab coordinate system is ...

z
: vertical direction

x
: horizontal, parallel to the projection of the beam center at the
sample position on the horizontal plane
: roghly pointing from *B&ouml;ttstein* towards *Villigen*

y
: horizontal, forming a right-handed system with x and z 
: roughly pointing from *Berg* towards *Aare*

The special situation with an incoming beam on the sample with
potentially a divergence of up to 1.5^o^ requires a differentiation 
between sample tilt ($\omega$) and angle of incidence ($\alpha_i$).
In the case of specular reflectivity we define $\theta := \alpha_i =
\alpha_f$. 

We define $\omega$ as the tilt angle of the sample with respect
to the horizon of the instrument. I.e. for liquid surfaces it is
allways zero!. 


*to be checked for consistency:* 
The position of the detector $\delta$ (realised by 2 translations in $x$ and 
$z$ direction), its tilt towards the sample $-\delta$ and the 
detection position of the neutron on / in the detector described by $\Delta\delta$
result in an angle $\delta + \Delta\delta = \omega + \theta$ relative to the
instrument horizon. Thus
$$ \theta = \delta + \Delta\delta - \omega$$

![Geometry and angles used for the focused beam set-up. The blue line
represents the horison, the black line the inclination $\omega$ of the sample surface
relative to the horizon. The gray bar represents the detector located
at an angle $\delta$. The red neutron trajectory hits the detector with an
offset $\Delta\delta$.](angles.eps)
 

## Set-up options

### polarisation

Amor is (will be) equipped with a spiral-shaped neutron polariser.
It is located before the first chopper within a vacuum housing in
the neutron guide bunker.

The polariser is mounted on a vertical translation system which 
allows for 3 settings: 

- **polarised**: the polariser with a frame overlap filter coating
  is in the beam
- **unpolarised**: the frame overlap filter is in the beam
- **open**: the neutron beam passes in between polariser and
  frame overlap filter

### liquid surfaces

(not yet installed)

The neutron beam can be inclined downward onto a horizontal (liquid)
surface by a mirror mounted on the optical bench right after the
first diaphragm.

The limited acceptance of this mirror means that only a small part of
the beam can be used. I.e. wide-divergence reflectometry is not
possible in this mode. The srrong off-specular scattering from a 
liquid surface anyhow asks for a low-divergent beam.

Within the divergence of the beam before the first diaphragm it is
possible to scan / vary the angle of incidence on the sample
without moving it vertically.
 
### beam divergence

The *natural* divergence of the beam at the sample position 
in the scattering plane is about 1.5 deg. It can be reduced by
closing the first slit behind the guide, D1. 

The independent movement of the blades allows to vary or scan 
the position of the slit opening across the divergent beam.
I.e. The angle of incidence on the sample can be changed without
tilting the sample.  

---

# cold start of the instrument

Short description of how to start hard- and software after the shut down 
or after a power failure.

## electronic racks

Starting the **motor control units** (MCU) in rack 1 and the **SPS** in rack 3 
has to be done by the repsective LIN staff. Presently that would be
Marccel Schild (MCU) and Roman Buerge (SPS).

Once the SPS is running, one has to set the gas flows for the detector and
probably for the flight tube:

```
SPS
 |- main
 |- Detector 
    |- gas flow 'Detector counting gas' to 7 sccm
    |- gas flow 'Flight tube Argon' to 10 sccm
```

## chopper

Starting the chopper is tricky, instable and time-consuming. Be patient.
Mind the order of the operations!

Situation: *device control unit* DCU (upper box, small display) and supervision 
(lower box with large screen) switched off

1. switch on DCU and boot  
   key switch to 'PC'
1. switch on supervision (back side), boot and wait for `AMOR GUI` to start
1. on the controll screen:

   1. klick on the top bar of *Chopper 1*
   1. tag *Chopper*: *DCU Command* to `Callibrate`  
      *accept-key*
   1. klick on the top bar of *Chopper 2*
   1. tag *Chopper*: *DCU Command* to `Callibrate`  
      *accept-key*
   1. wait for both choppers to report *positioning* or *ready*
      on the DCU screen
   1. klick on the top bar of *Chopper 1*
   1. tag *Chopper*: *DCU Command* to `Async Rotation`  
      *Speed* to 500 rpm  
      *accept-key*
   1. klick on the top bar of *Chopper 2*  
      *DCU Command*: `Sync. rotation`  
      *Phase*: -17  
      *Master*: `Master Chopper 1`  
      *Gear Ratio*: `1/1`  
      *accept-key*
   1. should be fine now!

## amor computer and instrument control

The instrument control computer `amor.psi.ch` (`172.28.65.60`, `PC14655`) 
is located in the cabin under the table on the right side.

1. switch on and boot
1. log in as user *amorlnsg* (24lns1)
1. start and control processes:  
   open a terminal window and switch user to *amor*:

   ```
   su - amor
   marche-gui
   ```
     
   in Marche gui:

   - klick `amor.psi.ch`
   - start *nicos*
   - start *amorIOC*
   - start *Histogram Memory*
   - start all *Kafka-to-Nexus...* processes 
   - start *Automatic File Sync*

1. start NICOS:
      
  - open a terminal window and enter  
    `nicos-gui &`
  - in the NICOS gui klick the tooth wheel  
    `Connect to server`  
    *User name*: `user`  
    *Password*: `24lns1`  
    *OK*

## detector 

see section *MBamor*

  
## data visualisation and reduction

as user `amorlnsg@amor.psi.ch` (i.e. when opening a new termminal on the amor computer)
perform the following steps:

```
source ~/eos/venv/bin/activate.csh
cd <working directory>
python ~/eos/events2histogram.py
display -update 2 e2h.png &
```

Now you can use `events2histogram.py` for fast visualisation and
`eos.py` for data reduction in THIS terminal window.

## sample environment

in NICOS gui:

```
|- Setup
|  |- Instrument
|  |  |- frappy       -> select
|  |  |- frappy_main  -> select
|  |- Apply
|
|- Instrument interaction
      |- output          -> type frappy('<device>')
```

e.g. potentiostat:  
`frappy('smamor')`
generates `frappy-main` with devices `se_smi` and `se_smv`

in terminal window:

```
sea
```

---

# Marche - daemon control

The daemons for NICOS, EPICS, `filewriter` and so on are controlled via
`Marche`. 

1. `marche` gui

> 1. log in as user `amor@amor.psi.ch`
> 1. call `marche-gui`

1. `marche` script ?



# the detector **MBamor**

## parameters

size: 14 blades, each 32 channels high and 64 channels wide  
sample detector distance (to blade tips): 4000 mm  
(distance blade tip to housing front: 33 mm  
inclination of the blades: 5.1 deg  
separation of blades: 10.5 mm or 0.15 deg  


## launching of electronics and hardware

1.  start server **det-efu02** in rack 1  
    (a monitor might be needed during booting to enable 'vnc')
1.  switch on **master module** in rack 1
1.  switch on **assistor**s on top of the detector tower  
    (probably just reconnect the LV cables)  
    check **Master/ring status**: should be 'not configured'
1.  ramp up **HV**

    | channel | name | U / V | I / $\mu$A |
    | ---     | ---  | ---:  | ---:       |
    | 0       | K1   | 1230  | 73         |
    | 1       | R1   | 1050  | 41         |
    | 2       | K2   | 1230  | 73         |
    | 3       | R2   | 1050  | 41         |
 
    (start from *off* position, not from *kill*)

## start processes on server **det-efu02**


The actions below can either be performed via `ssh` or `vnc`. In the latter case
(more challenging for the server graphic card) one might have to start the
**vncserver** first on the det-efu02:
    
> presently, vnc does not work - most likely because of RH8.
> use `ssh` instead!
>
> 1.  start `vncserver` (optional, should be launched at booting)
>
>     ```
>     cd essproj
>     ./startvnc
>     cd
>     ```
> 1.  change graphic settings for vnc
> 
>     ```
>     cd detg_git
>     ./display_vnc
>     cd
>     ```
> connecting via `vnc`
>
> -   via web browser from `amor` or `amor-dr`:
>     enter 'vnc://det-efu02:5901' in address bar
> -   via `boxes` on console  
>     `vnc://172.28.65.80:5901` mit passwd `essdaq`

1.  log into `det_efu02`:
    (e.g. from a terminal window on the computer `amor.psi.ch` as user `amor` or `amorlnsg`)

    ```
    ssh essdaq@172.28.65.80
    ``` 

    (password `essdaq`)

1.  start **ring controll**`

    ```
    cd slowctrl_rmm
    ./run_all.sh
    cd
    ```

    (the `poll` numbers should not exceed 16)

    the **master/ring status** shown on the master module should be green for rings 0, 1 and 2

1.  start slow controll pannel

    ```
    ./vmmdcs &
    ```

    `VMM3a` opens

    initialise slow control in `VMM3a` window

    1.  load configuration

        - *open folder* (lower left region)  
        - select 'Amordetector....' and *open*
        - *Load*  

        *FEN00\** shows up

    1.  initialise and start acquisition
 
        - *open communication* 
        - array on the left should show '4' for all hybrids

          > if not all hybrids show '4':
          > 
          > - *ACQ off*
          > - select *FEN0/1/2*
          > - *Warm init*
          > - *Send*
          > - *ACQ on*
          > 
          > escalation:
          > 
          > 1. *Hard reset* (lower right corner)
          > 1. power cycle LV (NOT front end assistors)

        - *send*
        - *ACQ on*
        - array on the left should show '5' for all hybrids


1.  start **event formation unit**s

    ```
    cd Desktop
    ./runall.bash
    cd
    ```

## start processes on `amor-dr`

### daqlite

start `daqlite` for *real*time images of detector output

```
cd Desktop
./daqlite.bash
cd
```

## start processes on `amor`
 
### graphana

to monitor the efu activity

1.  on amor:

    - open Firefox tab  
    - enter `https://172.28.65.80:3000` (or bookmark)  
    - select 'Freia'

2.  on `amor`:

    - open firefox browser

      address: `vnc://det-efu02:5901`  
      passwd: essdaq

    - open browser

      select AMORB2


## shut down

1.  **slow controll**: *ACQ off*
2.  switch off **front end assistor**s
3.  in random order:

    -   close slow controll
    -   close daqlite
    -   stop **efu** on det-efu02:

        ```
        cd Desktop
        ./stopall.bash
        cd
       ```

4.  shut down det-efu02  
    (can be powerd off without shutting down)
        
## configuration options

### chopper signal

The default setting is that the chopper trigger signal is used as a timing signal.
I.e. without chopper TTL signal, no data acqusition.  
This can be changed on 'det-efu02' in the script `.../run_all.sh` at approximately line
30:  
toggle the comments on the entries  
`...`  
and   
`...`  

## maintanance

### temperature

check every couple of weeks:  
*VMM3a / I2C / measure*

### uptime

The **VMM**s have a limited life time (approx 5 years). It is thus recommended to
power them down if not used for more than about a week:
*VMM3a / ACQ off*
disconnect LV power cable 

### trouble shooting

-   **dark lines
    *VMM3a / ACQ off*
    *VMM3a / FEN\* / Warm init* 
    *VMM3a / Send*
    *VMM3a / ACQ on*

-   **dark area**
    on detector image but all hybrids show '5' (green):  
    restart rings: 
    `det-efu02> .../run_all.sh`

-   **no data**

    -   chopper TTL signal might be missing
    -   *VMM3a / ... * all on '5' (green):
        restart rings:  
        `det-efu02> .../run_all.sh`
 
# trouble shooting

## (re)start services

- log in as `amor@amor.psi.ch`:
- open `marche-gui` (unless already open)
- check for status of

  - EPICS 
  - NICOS 
  - `chopper.service`
  - `histogrammer.service`
  - `kafka-to-nexus.target`

  and probably *stop* / *restart* / *start* the service(s)

- I do not remember what this is for:

  ```
  telnet localhost 20000
  ^x
  ^[
  quit
  ```

# kafka stuff

## check kafka input from `efu`

on `amor`:

- `source /home/software/virtualenvs/kafka-tools/bin/activate.csh`
- `kafka-spy -b linkafka01:9092 -t amor_ev44 -C`

# NICOS commands

## start / stop NICOS

in case the Selene pitches MCUs are active, one can deactivate these with
 
    amor> caput SQ:AMOR:SEL2.AOUT M714=0
    amor> caput SQ:AMOR:SEL2.AOUT M814=0

The **gui** can be launched also by the `amorlnsg` from the console:
and to launch the gui

    /home/amor/bin/nicos-gui

with

    user: user 
    password: 24lns1

## devices and movements

In NICOS everything is a device! 

`ListDevices()`

: shows all available devices 

`read(<device>)`

: reads the value of a *device* 

`maw(<device>, <value> [, <device2>, <value2> ...])`

: (move and wait) moves a *device* to *value* 

`move(<device>, <value>)`

: starts moving a *device* to *value* without waiting 

`wait(<dev1>, <dev2>, ...)`

: waits for the listed *devices* to finish 

`status(<device>)`

: queries the status of a *device*

`stop()`

: stops devices, also Ctr-C, Ctrl-C

NICOS devices can have parameters 

`ListParams(<device>)`

: lists the accessible parameters of a *device* 

`<device>.<parameter>`

: queries the value of a *parameter* of *device*

`<device>.<parameter>=<value>`

: sets a new *value* for a *parameter* of a *device*

`ListMethods(<device>)`

: lists the methods of a *device*

## measurement

`count(<mode>=<preset>)`

: count with *mode* 

  - `t` = time for *preset* seconds
  - `m` = monitor for *preset* monitor counts 

`scan(<device>, <start>, <step>, <np>, <mode>=<preset>)`

: scan the *device* from *start* with *np* steps of width *step* with the
counting time defined by *mode* / *preset*  

`scan(<device>, [<p1>, <p2>, <p3>, <p4>, ... ],  <mode>=<preset>)`

: scans points *p1, p2, ...* as given in the list 

`cscan(dev, center, step, np, t=2)`

: center-scan: scans *device* around a *center* with *np* steps of width *step* to
either side.

## batch files

Batch files are written in python.
You can use all of python plus the NICOS command set 

`run(<scriptname>)`

: *run* the script *scriptname*, which is either relative to 
`Exp.scriptpath` or has an absolute path 

`sim(<scriptname>)`

: *sim*ulate the script *scriptname*

## raw data

Location of the raw data files:

    /home/amor/data/<year>/<proposalname>

Location of lof files:

    /home/amor/nicos/log

