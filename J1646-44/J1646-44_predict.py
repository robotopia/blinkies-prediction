import numpy as np
import astropy.units as u
from astropy.time import Time, TimezoneInfo

#-------------------------------------------------#
# Input parameters and settings (user may change) #
#-------------------------------------------------#

# Change the range of desired MJD prediction here
MJD_start = Time(60250, format='mjd')
MJD_stop = Time(60260, format='mjd')

# Set the reference epoch (MJD)
T0 = Time(59391.3567, format='mjd')

# Output timezone
output_timezone = TimezoneInfo(utc_offset=5.5*u.hour)

# Offset reported times by some custom amount (e.g. half an hour before the ingress/egress times)?
# Set to None or 0.0 for no offset
offset = -0.5 * u.hour

# Pre-determined orbital properties
Pb = 5.267104 * u.hour
ingress_epoch = T0 + 0.0502*Pb
egress_epoch = T0 + 0.4497*Pb

#---------------------------------------#
# Calculations (user should not change) #
#---------------------------------------#

# Calculate the number of orbital periods since the assumed egress/ingress epochs that fall within the given MJD range
egress_rotation_start = int(np.ceil((MJD_start - egress_epoch)/Pb))
egress_rotation_stop = int(np.ceil((MJD_stop - egress_epoch)/Pb))
ingress_rotation_start = int(np.ceil((MJD_start - ingress_epoch)/Pb))
ingress_rotation_stop = int(np.ceil((MJD_stop - ingress_epoch)/Pb))

# Calculate the actual ingress/egress times
egresses = np.arange(egress_rotation_start, egress_rotation_stop)*Pb + egress_epoch
ingresses = np.arange(ingress_rotation_start, ingress_rotation_stop)*Pb + ingress_epoch

# Apply the offset
if offset is not None:
    egresses += offset
    ingresses += offset

# Convert to specified timezone and print out
print(f"Egresses {offset.to('hour'):+f}:")
for egress in egresses:
    print(egress.to_datetime(timezone=output_timezone))

print(f"Ingresses {offset.to('hour'):+f}:")
for ingress in ingresses:
    print(ingress.to_datetime(timezone=output_timezone))

