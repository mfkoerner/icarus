# This is the default configuration file. Copy it to config.py and set your specific configuration there.


# SET YOUR MATERIALS PROJECT API KEY HERE
matprojapi = "YOUR_pymatgen_API_KEY_HERE"


# SET A PATH TO YOUR MATPLOTLIB STYLESHEET HERE
plot_style = "path_to_a_matplotlib_stylesheet"


# CHANGE THIS FOR YOUR PARTICULAR SHARED DROPBOX LOCATION
shared_loc = 'YOUR_SHARED_LOCATION'


# CHANGE THIS FOR YOUR PARTICULAR LOCAL ARCHIVE LOCATION
local_archive = 'YOUR_LOCAL_ARCHIVE_LOCATION'


# Default list of vasp directories for automated runs
VASP_directories = ['static', 'band', 'pbesoc', 'SOCband', 'HSESOC', 'SHband', 'ktest', 'entest', 'absorb', 'DOS']


# Default path to status.txt for vasp runs (contact mitchell for explanation)
statuspath = 'path to status here'

# Install Location here
install_loc = 'install_location'

# mpi-run settings for vasp
VASP_nodes = 1
VASP_ppn   = 16
VASP_walltime = 1
VASP_symprec = 1e-5
