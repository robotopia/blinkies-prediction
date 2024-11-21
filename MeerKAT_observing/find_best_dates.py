import numpy as np
import astropy.units as u
from astropy.time import Time, TimezoneInfo
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.time import Time

from astroplan import Observer
from astroplan import AltitudeConstraint, PhaseConstraint
from astroplan import is_observable, is_always_observable, months_observable
from astroplan import FixedTarget
from astroplan.plots import plot_altitude
from astroplan import PeriodicEvent

from pytz import timezone
import matplotlib.pyplot as plt

mkt = Observer.at_site("salt")

# We received 1 hour of observing time per source
# MeerKAT documentation says up to 30-minute integrations are possible
# Pulsar-searching rule-of-thumb is to not search over more than 0.1*Pb, i.e. a max of 30 minutes for our P=5-hour sources
# But splitting the blocks more might gain us more chances, so good options for this value are 0.25 (4x15-min observations per source) and 0.5 (2x30-min observations per source)
blocksize = 0.5 # hours
nobs = int(1/blocksize)
# Also print the 'ok' days not just the ideal days (currently disabled)
#allowTolerable = False

coords = {"J1646": FixedTarget(SkyCoord("16 46 22.7", "-44 05 41", frame="fk5", unit=(u.hour, u.deg)), name="J1646"),
          "J1723" : FixedTarget(SkyCoord("17 23 38.5", "-33 42 00", frame="fk5", unit=(u.hour, u.deg)), name="J1723"),
          "J1728" : FixedTarget(SkyCoord("17 28 12.12", "-46 08 01.49", frame="fk5", unit=(u.hour, u.deg)), name="J1728"),
          "J1734" : FixedTarget(SkyCoord("17 34 33.9", "-28 07 07", frame="fk5", unit=(u.hour, u.deg)), name="J1734"),
          "J1740" : FixedTarget(SkyCoord("17 40 16.16", "-26 50 28.77", frame="fk5", unit=(u.hour, u.deg)), name="J1740"),
          "J1752" : FixedTarget(SkyCoord("17 52 26.2", "-30 33 43", frame="fk5", unit=(u.hour, u.deg)), name="J1752"),
          "J1815" : FixedTarget(SkyCoord("18 15 56.8", "-14 16 34", frame="fk5", unit=(u.hour, u.deg)), name="J1815")}

# Orbital periods
Pbs = {"J1646" : 5.26703*u.hour,
       "J1723" : 22.68495*u.hour,
       "J1728" : 5.04983*u.hour,
       "J1734" : 10.92946*u.hour,
       "J1740" : 16.50614*u.hour,
       "J1752" : 17.85847*u.hour,
       "J1815" : 9.81969*u.hour}

# T0 = middle of eclipse, in MJD
# I have scaled J1646's T0 to add 0.26 of orbital phase based on Fig 3
# I have pulled back J1740's T0 to 0.9 of orbital phase since it has a very slow egress
T0s = {"J1646" : Time(59391.3567+0.26*(Pbs["J1646"].to(u.day).value), format='mjd', scale='utc', location=mkt.location),
       "J1723" : Time(59789.566010035276, format='mjd', scale='utc', location=mkt.location),
       "J1728" : Time(60055.87090177478, format='mjd', scale='utc', location=mkt.location),
       "J1734" : Time(59796.055897211045, format='mjd', scale='utc', location=mkt.location),
       "J1740" : Time(60055.87090177478-0.1*(Pbs["J1740"].to(u.day).value), format='mjd', scale='utc', location=mkt.location),
       "J1752" : Time(59807.08632836161, format='mjd', scale='utc', location=mkt.location),
       "J1815" : Time(59740.04800876778, format='mjd', scale='utc', location=mkt.location)}

# duty cycle of ON phase
duty_cycles = {"J1646" : 0.5,
       "J1723" : 0.7,
       "J1728" : 0.6,
       "J1734" : 0.7,
       "J1740" : 0.2,
       "J1752" : 0.4,
       "J1815" : 0.5}


# Astroplan Periodic Event objects
pe = {}
for obj in T0s:
    pe[obj] = PeriodicEvent(epoch=T0s[obj], period=Pbs[obj])

# Altitude constraints at MeerKAT -- 15 deg is in the OPT, but I'm getting conflicting results from astroplan, so make it higher
ac = AltitudeConstraint(min=30*u.deg)

# Observing constraints (right part of orbital phase, and above the horizon)
# Note that astroplan's default is quite loose constraints -- if your observation has ANY 'on' time, then it counts as OK
# E.g. you could start observing JUST before egress and then most of the observation would be useless
# So to make sure the WHOLE observation fits within the orbital phase constraint, you need to shrink the window by the blocksize
# Which to turn into a phase constraint, is blocksize/Pb

time_window = {}
constraints = {}
print("Source  Min Orb Phase  Max Orb Phase")
for obj in T0s:
    time_window[obj] = Pbs[obj]*duty_cycles[obj] - blocksize*u.hour
    constraints[obj] = [PhaseConstraint(pe[obj], min=0.5 - (duty_cycles[obj]/2) + blocksize*u.hour/Pbs[obj], max=0.5 + (duty_cycles[obj]/2) - blocksize*u.hour/Pbs[obj]), ac]
    print(obj, .5 - (duty_cycles[obj]/2) + blocksize*u.hour/Pbs[obj], 0.5 + (duty_cycles[obj]/2) - blocksize*u.hour/Pbs[obj])

# ordering: by the time window of 'on'-ness
order = {k: v for k, v in sorted(time_window.items(), key=lambda item: item[1])}


#-------------------------------------------------#
# Input parameters and settings (user may change) #
#-------------------------------------------------#


# Change the range of desired MJD prediction here
MJD_start = Time("2024-12-31T00:00:00", format='isot', scale="utc", location=mkt.location)
#MJD_stop = Time("2025-11-15T00:00:00", format='isot', scale="utc", location=mkt.location)
# For debugging -- and in the hopes they can schedule quickly
MJD_stop = Time("2025-01-31T00:00:00", format='isot', scale="utc", location=mkt.location)

#---------------------------------------#
# Calculations (user should not change) #
#---------------------------------------#

# For each day, calculate whether all of the sources are observable or if there is a showstopper:
ok = []
for mjd in np.arange(MJD_start.mjd, MJD_stop.mjd):
    keep = True
    for obj in order:
        if not is_observable(constraints[obj], mkt, coords[obj], time_range=[Time(mjd, format='mjd', scale='utc'), Time(mjd+1, format='mjd', scale='utc')]):
            print(f"Can't observe {obj} on {mjd}")
            keep = False
    if keep:
        print(f"All sources available on {mjd}")
        ok.append(mjd) 

# Find out when the start and stop times are for a given day

obslength = []
for mjd in ok:
    t = Time(mjd, format='mjd', scale='utc')
    start = False
    stop = False
    for hour in range(0, 24):
        observability = [is_observable(constraints[obj], mkt, coords[obj], time_range=[Time(mjd+(hour/24), format='mjd', scale='utc'), Time(mjd+((hour+1)/24.), format='mjd', scale='utc')])[0] for obj in order]
        if not start and np.any(observability):
             start = hour
        if start and not stop and not np.any(observability):
             stop = hour
    print(f"For {t.isot}, observability starts at {start} and ends at {stop}")
    if stop - start > 8:
    
# Create a detailed schedule by observing each source for blocksize minutes and moving on to the next source with the fewest observations so far
        scheduled_times = []
        scheduled_objects = []
        #minobs = 0
# Number of (blocksize min) observations taken so far (zero for all objects)
        observed = {"J1646" : 0,
               "J1723" : 0,
               "J1728" : 0,
               "J1734" : 0,
               "J1740" : 0,
               "J1752" : 0,
               "J1815" : 0}

        skipsource = None
        for hour in np.arange(start, stop, blocksize):
             blockAvailable = True
             for obj in order:
                  if is_observable(constraints[obj], mkt, coords[obj], time_range=[Time(mjd+(hour/24), format='mjd', scale='utc'), Time(mjd+((hour+blocksize)/24.), format='mjd', scale='utc')]) and observed[obj] < nobs and blockAvailable and obj != skipsource:
                  #if is_observable(constraints[obj], mkt, coords[obj], time_range=[Time(mjd+(hour/24), format='mjd', scale='utc'), Time(mjd+((hour+blocksize)/24.), format='mjd', scale='utc')]) and observed[obj] < nobs and observed[obj] <= minobs and blockAvailable:
#                       print(obj, Time(mjd+(hour/24), format='mjd', scale='utc').isot, pe[obj].phase(Time(mjd+(hour/24), format='mjd', scale='utc')), pe[obj].phase(Time(mjd+((hour+blocksize)/24.), format='mjd', scale='utc')))
                       scheduled_times.append(Time(mjd+(hour/24), format='mjd', scale='utc'))
                       scheduled_objects.append(obj)
                       observed[obj] += 1
                       #minobs = min(observed.values())
                       blockAvailable = False
                       skipsource = obj

        if len(scheduled_times) == nobs*len(T0s):
            print(f"{t.isot} is an ideal day, with a schedule of {len(scheduled_times)} {blocksize*60}-minute observing blocks that starts at {start}h, ends at {stop}h, and iterates through the sources in the following order:")
            print("Time/UTC source orb_phase")
            phases = []
            for t,source in zip(scheduled_times, scheduled_objects):
                print(t.isot, source, pe[source].phase(t))
                phases.append(pe[source].phase(t))
            with open("schedule_options.txt", "a") as f:
                f.write(f"{scheduled_times[0].isot} {np.nanmean(phases)} {np.nanstd(phases)}\n")
#        elif len(schedule) > nobs*len(T0s) - 4 and allowTolerable:
#            print(f"{t.isot} is an tolerable day, with a schedule of {len(schedule)} {blocksize*60}-minute observing blocks that starts at {start}h, ends at {stop}h, and iterates through the sources in the following order:")
#            print(schedule)
        else:
            print(f"Couldn't spread the sources out cleanly across {t.isot}; schedule was only {len(scheduled_times)*blocksize} hours long")
    else:
        print(f"Only {stop - start} hours of useful time available on {t.isot}.")

