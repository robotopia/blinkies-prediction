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

coords = {"J1646": FixedTarget(SkyCoord("16 46 22.7", "-44 05 41", frame="fk5", unit=(u.hour, u.deg)), name="J1646"),
          "J1723" : FixedTarget(SkyCoord("17 23 38.5", "-33 42 00", frame="fk5", unit=(u.hour, u.deg)), name="J1723"),
          "J1728" : FixedTarget(SkyCoord("17 28 12.12", "-46 08 01.49", frame="fk5", unit=(u.hour, u.deg)), name="J1728"),
          "J1734" : FixedTarget(SkyCoord("17 34 33.9", "-28 07 07", frame="fk5", unit=(u.hour, u.deg)), name="J1734"),
          "J1740" : FixedTarget(SkyCoord("17 40 16.16", "-26 50 28.77", frame="fk5", unit=(u.hour, u.deg)), name="J1740"),
          "J1752" : FixedTarget(SkyCoord("17 52 26.2", "-30 33 43", frame="fk5", unit=(u.hour, u.deg)), name="J1752"),
          "J1815" : FixedTarget(SkyCoord("18 15 56.8", "-14 16 34", frame="fk5", unit=(u.hour, u.deg)), name="J1815")}

# T0 = middle of eclipse, in MJD
# I have scaled J1646's T0 to 0.26 of orbital phase based on Fig 3
T0s = {"J1646" : Time(59391.3567+0.26*(5.26703/24.), format='mjd', scale='utc', location=mkt.location),
       "J1723" : Time(59789.566010035276, format='mjd', scale='utc', location=mkt.location),
       "J1728" : Time(60055.87090177478, format='mjd', scale='utc', location=mkt.location),
       "J1734" : Time(59796.055897211045, format='mjd', scale='utc', location=mkt.location),
       "J1740" : Time(60056.868761778795, format='mjd', scale='utc', location=mkt.location),
       "J1752" : Time(59807.08632836161, format='mjd', scale='utc', location=mkt.location),
       "J1815" : Time(59740.04800876778, format='mjd', scale='utc', location=mkt.location)}

# Orbital periods
Pbs = {"J1646" : 5.26703*u.hour,
       "J1723" : 22.68495*u.hour,
       "J1728" : 5.04983*u.hour,
       "J1734" : 10.92946*u.hour,
       "J1740" : 16.50614*u.hour,
       "J1752" : 17.85847*u.hour,
       "J1815" : 9.81969*u.hour}

# duty cycle of ON phase
duty_cycles = {"J1646" : 0.5,
       "J1723" : 0.8,
       "J1728" : 0.6,
       "J1734" : 0.7,
       "J1740" : 0.2,
       "J1752" : 0.4,
       "J1815" : 0.5}

# Astroplan Periodic Event objects
pe = {}
for obj in T0s:
    pe[obj] = PeriodicEvent(epoch=T0s[obj], period=Pbs[obj])

# Altitude constraints at MeerKAT
ac = AltitudeConstraint(min=15*u.deg)

# Observing constraints (right part of orbital phase, and above the horizon)
constraints = {}
for obj in T0s:
    constraints[obj] = [PhaseConstraint(pe[obj], min=0.5 - (duty_cycles[obj]/2), max=0.5 + (duty_cycles[obj]/2)), ac]

#-------------------------------------------------#
# Input parameters and settings (user may change) #
#-------------------------------------------------#


# Change the range of desired MJD prediction here
MJD_start = Time("2024-11-30T00:00:00", format='isot', scale="utc", location=mkt.location)
MJD_stop = Time("2025-11-15T00:00:00", format='isot', scale="utc", location=mkt.location)
# For debugging
#MJD_stop = Time("2024-12-15T00:00:00", format='isot', scale="utc", location=mkt.location)


#---------------------------------------#
# Calculations (user should not change) #
#---------------------------------------#

# For each day, calculate whether all of the sources are observable or if there is a showstopper:
ok = []
for mjd in np.arange(MJD_start.mjd, MJD_stop.mjd):
    keep = True
    for obj in T0s:
        if not is_observable(constraints[obj], mkt, coords[obj], time_range=[Time(mjd, format='mjd', scale='utc'), Time(mjd+1, format='mjd', scale='utc')]):
            print(f"Can't observe {obj} on {mjd}")
            keep = False
    if keep:
        ok.append(mjd) 

# Find out when the start and stop times are for a given day

obslength = []
for mjd in ok:
    start = False
    stop = False
    for hour in range(0, 24):
        observability = [is_observable(constraints[obj], mkt, coords[obj], time_range=[Time(mjd+(hour/24), format='mjd', scale='utc'), Time(mjd+((hour+1)/24.), format='mjd', scale='utc')])[0] for obj in T0s]
        if not start and np.any(observability):
             start = hour
        if start and not stop and not np.any(observability):
             stop = hour
    print(f"For {mjd}, observability starts at {start} and ends at {stop}")
    if stop - start > 8:
    
# Create a detailed schedule by observing each source for 15 minutes and moving on to the next source with the fewest observations so far
        schedule = []
        maxobs = 1
        minobs = 0
# Number of (15-min) observations taken so far (zero for all objects)
        observed = {"J1646" : 0,
               "J1723" : 0,
               "J1728" : 0,
               "J1734" : 0,
               "J1740" : 0,
               "J1752" : 0,
               "J1815" : 0}

        for hour in np.arange(start, stop, 0.25):
             for obj in T0s:
                  if is_observable(constraints[obj], mkt, coords[obj], time_range=[Time(mjd+(hour/24), format='mjd', scale='utc'), Time(mjd+((hour+0.25)/24.), format='mjd', scale='utc')]) and observed[obj] <= 4 and observed[obj] <= minobs:
                       schedule.append(obj)
                       observed[obj] += 1
                       minobs = min(observed.values())

        if len(schedule) > 25:
            print(f"{mjd} is an ideal day, with a schedule that starts at {start}, ends at {stop}, and iterates through the sources in the following order:")
            print(schedule)
