Water for Malkmus Mo
--------------------

The fields output for this profile are listed in order below with format strings, units and description the following information: :

global_iso_id
-------------
Data type: int
Units: [dimensionless]
Description: Unique integer ID of a particular isotopologue: every global isotopologue ID is unique to a particular species, even between different molecules. The number itself is, however arbitrary.

nu
--
Data type: float
Units: cm-1
Description: Transition wavenumber

sw
--
Data type: float
Units: cm-1/(molec.cm-2)
Description: Line intensity, multiplied by isotopologue abundance, at T = 296 K

gamma_self
----------
Data type: float
Units: cm-1.atm-1
Description: Self-broadened HWHM at 1 atm pressure and 296 K

n_air
-----
Data type: float
Units: [dimensionless]
Description: Temperature exponent for the air-broadened HWHM

elower
------
Data type: float
Units: cm-1
Description: Lower-state energy
