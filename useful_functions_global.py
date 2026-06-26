import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings
import seaborn as sns

def flux2mag(flux):
    zeropoint = 31.4
    index = flux.index if hasattr(flux, 'index') else None
    flux = np.asarray(flux, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        mag = -2.5 * np.log10(flux) + zeropoint
    mag[~np.isfinite(mag)] = np.nan
    if index is not None:
        return pd.Series(mag, index=index)
    return mag

def fluxerr2magerr(flux, flux_err):
    flux = np.ma.filled(np.ma.asarray(flux, dtype=float), fill_value=np.nan)
    flux_err = np.ma.filled(np.ma.asarray(flux_err, dtype=float), fill_value=np.nan)
    with np.errstate(invalid='ignore', divide='ignore'):
        magerr = (2.5 / np.log(10)) * (flux_err / flux)
    magerr[~np.isfinite(magerr)] = np.nan
    return magerr


def get_tract(field):
    #for DP1, from https://portal.nersc.gov/cfs/lsst/dp1/contributed-notebooks/DP1_Detector_Visits_NB1.html
    tract_dict={'47 Tuc' : [531, 532, 453, 454],
            'ECDFS' : [4848, 4849, 5063, 4636, 4637, 4638, 4847, 4850, 5061, 5062, 5064,
       5065, 5279, 5280, 5281, 5282], #Euclid's "Euclid Deep Field Fornax" encompasses Extended Chandra Deep Field South
            'EDFS' : [2078, 2079, 2080, 2232, 2233, 2234, 2235, 2236, 2237, 2392, 2393,
       2394, 2395, 2396, 2397, 2557, 2558, 2559, 2560, 2561, 2562, 2728,
       2729, 2730, 2731],
            'Fornax' : [4016, 4017, 4218, 4217],
            'FDSG' : [4016, 4017, 4218, 4217], #same as Fornax idk how it'll be referred to
            'Rubin_SV_095-25' : [5525, 5526],
            'Rubin_SV_38_7' : [10463, 10464, 10704],
            'LELF': [10464, 10221, 10222, 10704, 10705, 10463], # Low Ecliptic Latitude Field / Rubin_SV_38_7
            'Seagull' : [7850, 7849, 7610, 7611],
            -1 : 'Not found'
       }
    key_return = -1
    for key in tract_dict.keys():
        if key == field:
            key_return = key
        else:
            continue
    return tract_dict[key_return]

def get_field(tract):
    field_dict={531 : '47 Tuc', 532 : '47 Tuc', 453 : '47 Tuc', 454 : '47 Tuc',
                5062 : 'ECDFS', 5063 : 'ECDFS', 5064 : 'ECDFS', 4848 : 'ECDFS', 4849 : 'ECDFS',
                2393 : 'EDFS', 2234 : 'EDFS', 2235 : 'EDFS', 2394 : 'EDFS',
                4016 : 'Fornax', 4017 : 'Fornax', 4218 : 'Fornax', 4217 : 'Fornax', #same as FDSG idk how it'll be referred to
                5525 : 'Rubin_SV_095-25', 5526 : 'Rubin_SV_095-25',
                10464 : 'LELF', 10221 : 'LELF', 10222 : 'LELF', 10704 : 'LELF', 10705 : 'LELF', 10463 : 'LELF', #same as Rubin_SV_38_7 I think
                7850 : 'Seagull', 7849 : 'Seagull', 7610 : 'Seagull', 7611 : 'Seagull',
            -1 : 'Not found'
       }
    key_return = -1
    for key in field_dict.keys():
        if key == tract:
            key_return = key
        else:
            continue
    return field_dict[key_return]


## PLOTS

params = {'legend.fontsize': 'x-large',
          'figure.figsize': (15, 5),
          'axes.labelsize': 'x-large',
          'axes.titlesize':'x-large',
          'xtick.labelsize':'x-large',
          'ytick.labelsize':'x-large'}
pad=15
size=1

def color_magnitude(star_df, band1, band2, colors, color_label, selection_label, title, presentation_mode = False, x_lim = (-1, 4), y_lim = (30, 18), ax = None, save = False, plots_path = None, filename = None):
    """
    Plots color-magnitude diagram, band1 vs band1-band2

    Parameters
    ----------
    star_df : Astropy Table
        Stellar catalog
    band1 : string
        Photometry band to be plotted on y-axis
    band2 : string
        Photometry band for the x-axis
        function wil plot {band}_psfFlux
    colors : Table column
        Colorbar data
    color_label : string
        Colorbar label
    title : string
        Title of plot/subplot
    ax (optional) : default None, else plt axes object
        Axes to plot subplots on
        If none, will set subplot axes (1,1)
    save (optional) : default False
        File names have form 'colormag_{title}.png'
    plots_path (optional) : default None, else string
        Where to save

    Returns
    -------
    Pretty plot
    """

    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(7,5))
        ax = axes

    band1_mag = flux2mag(star_df[f'{band1}_psfFlux'])
    band2_mag = flux2mag(star_df[f'{band2}_psfFlux'])
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='.*colormapping.*')
        warnings.filterwarnings("ignore", message='.*labels.*')
        if presentation_mode == True:
            ax.hist2d(band1_mag - band2_mag, band1_mag, bins=200,
                       range = [[x_lim[0],x_lim[1]],[y_lim[1], y_lim[0]]],
                       cmin=1, cmap = colors, label = selection_label)
        else:
            _ = ax.scatter(band1_mag - band2_mag, band1_mag,
                    c = colors, vmin=0, vmax=1,
                    label = selection_label,
                    s = size, cmap = 'jet_r')
    if selection_label[0] != '_':
        ax.legend()
    if color_label is not None:
        plt.colorbar(_,label=color_label)
    ax.set_title(title, pad=pad)
    ax.set(xlabel = f"{band1} - {band2}", ylabel = f"{band1}")
    ax.set(xlim = x_lim, ylim = y_lim)
    plt.tight_layout()
    if save == True:
        if filename is None:
            title = title.replace(' ', '').replace('-', '_').lower()
            filename = title
        plt.savefig(plots_path + f'/colormag/{filename}.png')


def color_magnitude2(star_df, pltL_dict, pltR_dict, title = None, save = False, plots_path = None, filename = None):
    """
    Plots 2 subplots of color-magnitude diagrams, band1 vs band1-band2

    Parameters
    ----------
    star_df : Astropy Table
        Stellar catalog
    pltL/R_dict : dictionary
        Params for left/right subplot
        {band1 : , band2 : , colors : , color_label : , title : }
        See color_magnitude for definitions of above
    title (optional) : default None, else string,
        Super title for plot
    save (optional) : default False
        File names have form 'colormag_{title}.png'
    plots_path (optional) : default None, else string
        Where to save

    Returns
    -------
    Pretty plot
    """

    fig, axes = plt.subplots(1,2, figsize = (14,5))
    axes = axes.flatten()
    dict_list = [pltL_dict, pltR_dict]

    for i, ax in enumerate(axes):
        color_magnitude(star_df,
                        dict_list[i]['band1'], dict_list[i]['band2'],
                        dict_list[i]['colors'], dict_list[i]['color_label'], dict_list[i]['title'],
                        ax = ax)

    if title is not None:
        fig.suptitle(title)
    plt.tight_layout()
    if save == True:
        if filename is None:
            title = title.replace(' ', '').replace('-', '_').lower()
            filename = title
        plt.savefig(plots_path + f'/colormag/{filename}.png')


def color_color(band_list, colors, color_label, selection_label, title, presentation_mode = False, y_lim = None, x_lim = None, ax = None, save = False, plots_path = None, filename = None):
    """
    Plots color-color, r-i vs g-r for different magnitudes

    Parameters
    ----------
    band_list : list of tuples
        1st entry of tuple is string labeling band
        2nd entries are the band data
    colors : Table column
        Data to be used for colorbar
    color_label : string
        Colorbar label
    title : string
        Plot title
    ax (optional) : default None, else plt axes object
        Axes to plot subplots on
        If none, will set subplot axes (1,1)
    save (optional) : default False
        File names have form 'colorcolor_{title}.png'
    plots_path : default None, else string
        Where to save

    Returns
    -------
    Pretty plot
    """

    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(9,8))
        ax = axes

    band_x1, band_x2, band_y1, band_y2 = band_list
    x1_mag = band_x1[1]
    x2_mag = band_x2[1]
    y1_mag = band_y1[1]
    y2_mag = band_y2[1]
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='.*colormapping.*')
        warnings.filterwarnings("ignore", message='.*labels.*')
        if presentation_mode == True:
            ax.hist2d(x1_mag - x2_mag, y1_mag - y2_mag, bins=200,
                       range = [[x_lim[0],x_lim[1]],[y_lim[0], y_lim[1]]],
                       cmin=1, cmap = colors, label = selection_label)
        else:
            _ = ax.scatter(x1_mag - x2_mag, y1_mag - y2_mag,
                    c = colors, vmin=0, vmax=1,
                    label = selection_label,
                    s = size, cmap = "jet_r")
    if selection_label[0] != '_':
        ax.legend()
    if color_label is not None:
        plt.colorbar(_,label=color_label)
    ax.set_title(title, pad=pad)
    ax.set(ylabel = f'{band_y1[0]} - {band_y2[0]}',
           xlabel = f'{band_x1[0]} - {band_x2[0]}')
    if y_lim is not None:
        ax.set(ylim=y_lim)
    if x_lim is not None:
        ax.set(xlim=x_lim)
    plt.tight_layout()
    if save == True:
        if filename is None:
            title = title.replace(' ', '').replace('-', '_').lower()
            filename = title
        plt.savefig(plots_path + f'/colorcolor/{filename}.png')


def star_gal_sep(df, seperator, colors, color_label, selection_label, title,
                 presentation_mode = False,
                 ax = None, colorbar_limits = (0,1), y_bounds=(-0.1, 0.6), x_bounds = (18,28),
                 line_plt = None, save = False, plots_path = '', filename = None):
    '''
    Plots the difference between PSF and cModel flux across magnitudes
    color-coded by star-galaxy classifiers
    Shows magnitude where things start getting confused

    Parameters
    ----------
    df : astropy table
        Table of all sources
    survey_sep : string
        Survey used for the separation, determines the y-axis
    colors : table column, potentially a
        Data to be used for colorbar
    color_label : string
        Colorbar label
    title : string
        Plot title
    ax (optional) : default None, else plt axes object
        Axes to plot subplots on
        If none, will set subplot axes (1,1)
    y_bounds (optional) : tuple of floats
        Argument for y limits
        Found it necessary sometimes to zoom and enhance on LSST selector
    save (optional) : default False
        If True the file will be saved
    plots_path (optional) : default empty string
        Where to save

    Returns
    -------
    Pretty plot
    '''
    if ax == None:
        fig, axes = plt.subplots(1,1, figsize=(7,5))
        ax = axes

    i_mag = flux2mag(df['i_psfFlux'])
    if seperator == 'i psf - cmodel':
        i_mag_cmodel = flux2mag(df["i_cModelFlux"])
        y = i_mag - i_mag_cmodel
        y_label = f'LSST i_psfFlux - i_cModelFlux'
    elif seperator == 'i psf / cmodel':
        i_mag_cmodel = flux2mag(df["i_cModelFlux"])
        y = i_mag / i_mag_cmodel
        y_label = f'LSST i_psfFlux / i_cModelFlux'
    else:
        y = df[seperator]
        y_label = 'Euclid ' + seperator
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message='.*colormapping.*')
        warnings.filterwarnings("ignore", message='.*labels.*')
        if presentation_mode == True:
            if type(colors) == str:
                ax.hist2d(i_mag, y, bins=200,
                       range = [[x_bounds[0],x_bounds[1]],[y_bounds[0], y_bounds[1]]],
                       cmin=1, cmap = colors, label = selection_label)
            else:
                _ = ax.scatter(i_mag, y,
                       c = colors, vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                       label = selection_label,
                       s = size, cmap = "viridis_r")
                sns.kdeplot(x=i_mag, y=y, fill=False, color="k")
        else:
            _ = ax.scatter(i_mag, y,
                       c = colors, vmin=colorbar_limits[0], vmax=colorbar_limits[1],
                       label = selection_label,
                       s = size, cmap = "viridis_r")
    if selection_label[0] != '_':
        ax.legend()
    if line_plt is not None:
        ax.plot(line_plt[0],line_plt[1],line_plt[2])
    ax.set(xlabel = 'LSST i_psfFlux mag', ylabel = y_label)
    if presentation_mode == False:
        ax.set(ylim = y_bounds, xlim = x_bounds)
    ax.set_title(title, pad=pad)
    if color_label is not None:
        plt.colorbar(_,label=color_label)
    plt.tight_layout()

    if save == True:
        if filename is None:
            title.replace(' ', '').replace('-','_').lower()
            filename = title
        plt.savefig(plots_path + f'/stargalsep/{filename}.png')


#There is definitely a better way to import ugali tools, but can't install via pip here.
#My temporary workaround is to just copy and paste from their github https://github.com/DarkEnergySurvey/ugali
"""
Class for converting between sphere to image coordinates using map projections.

Based on Calabretta & Greisen 2002, A&A, 357, 1077-1122
http://adsabs.harvard.edu/abs/2002A%26A...395.1077C
"""
############################################################

class SphericalRotator:
    """
    Base class for rotating points on a sphere.

    The input is a fiducial point (deg) which becomes (0, 0) in rotated coordinates.
    """

    def __init__(self, lon_ref, lat_ref, zenithal=False):
        self.setReference(lon_ref, lat_ref, zenithal)

    def setReference(self, lon_ref, lat_ref, zenithal=False):

        if zenithal:
            phi = (np.pi / 2.) + np.radians(lon_ref)
            theta = (np.pi / 2.) - np.radians(lat_ref)
            psi = 0.
        if not zenithal:
            phi = (-np.pi / 2.) + np.radians(lon_ref)
            theta = np.radians(lat_ref)
            # psi = 90 corresponds to (0, 0)
            # psi = -90 corresponds to (180, 0)
            psi = np.radians(90.)

        cos_psi,sin_psi = np.cos(psi),np.sin(psi)
        cos_phi,sin_phi = np.cos(phi),np.sin(phi)
        cos_theta,sin_theta = np.cos(theta),np.sin(theta)

        self.rotation_matrix = np.array([
            [cos_psi * cos_phi - cos_theta * sin_phi * sin_psi,
             cos_psi * sin_phi + cos_theta * cos_phi * sin_psi,
             sin_psi * sin_theta],
            [-sin_psi * cos_phi - cos_theta * sin_phi * cos_psi,
             -sin_psi * sin_phi + cos_theta * cos_phi * cos_psi,
             cos_psi * sin_theta],
            [sin_theta * sin_phi,
             -sin_theta * cos_phi,
             cos_theta]
        ])
        
        self.inverted_rotation_matrix = np.linalg.inv(self.rotation_matrix)

    def cartesian(self,lon,lat):
        lon = np.radians(lon)
        lat = np.radians(lat) 
        
        x = np.cos(lat) * np.cos(lon)
        y = np.cos(lat) * np.sin(lon)
        z =  np.sin(lat)
        return np.array([x,y,z])
        

    def rotate(self, lon, lat, invert=False):
        vec = self.cartesian(lon,lat)

        if invert:
            vec_prime = np.dot(np.array(self.inverted_rotation_matrix), vec)
        else:        
            vec_prime = np.dot(np.array(self.rotation_matrix), vec)

        lon_prime = np.arctan2(vec_prime[1], vec_prime[0])
        lat_prime = np.arcsin(vec_prime[2])

        return (np.degrees(lon_prime) % 360.), np.degrees(lat_prime)

############################################################

######

# ADW: Unteseted dummy projection.
def cartesianSphereToImage(lon, lat):
    lon = lon - 360.*(lon>180)
    x,y = lon,lat
    return x,y

def cartesianImageToSphere(x,y):
    x = x - 360.*(x>180)
    lon,lat = x,y
    return lon,lat

############################################################

def angsep2(lon_1, lat_1, lon_2, lat_2):
    """
    Angular separation (deg) between two sky coordinates.
    """
    import healpy

    v10, v11, v12 = healpy.ang2vec(np.radians(90. - lat_1), np.radians(lon_1)).transpose()
    v20, v21, v22 = healpy.ang2vec(np.radians(90. - lat_2), np.radians(lon_2)).transpose()
    val = (v10 * v20) + (v11 * v21) + (v12 * v22)
    val = np.clip(val, -1., 1.)
    return np.degrees(np.arccos(val))

def angsep(lon1,lat1,lon2,lat2):
    """
    Angular separation (deg) between two sky coordinates.
    Borrowed from astropy (www.astropy.org)

    Notes
    -----
    The angular separation is calculated using the Vincenty formula [1],
    which is slighly more complex and computationally expensive than
    some alternatives, but is stable at at all distances, including the
    poles and antipodes.

    [1] http://en.wikipedia.org/wiki/Great-circle_distance
    """
    lon1,lat1 = np.radians([lon1,lat1])
    lon2,lat2 = np.radians([lon2,lat2])
    
    sdlon = np.sin(lon2 - lon1)
    cdlon = np.cos(lon2 - lon1)
    slat1 = np.sin(lat1)
    slat2 = np.sin(lat2)
    clat1 = np.cos(lat1)
    clat2 = np.cos(lat2)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    return np.degrees(np.arctan2(np.hypot(num1,num2), denominator))

############################################################

def gal2cel(lon, lat):
    """Convert from Galactic coordinates (deg) to coordinates in the
    Celestial Equatorial J2000 (deg) frame.
    
    Parameters:
    -----------
    lon : Galactic longitude (deg)
    lat : Galactic latitude (deg)

    Returns
    -------
    ra,dec : Right ascension and declination (deg,deg)

    """
    lat = np.radians(lat)
    sin_lat = np.sin(lat)
    cos_lat = np.cos(lat)

    lon = np.radians(lon)
    ra_gp = np.radians(192.85948)
    de_gp = np.radians(27.12825)
    lcp = np.radians(122.932)

    sin_lcp_lon = np.sin(lcp - lon)
    cos_lcp_lon = np.cos(lcp - lon)

    sin_d = (np.sin(de_gp) * sin_lat) \
            + (np.cos(de_gp) * cos_lat * cos_lcp_lon)
    ramragp = np.arctan2(cos_lat * sin_lcp_lon,
                            (np.cos(de_gp) * sin_lat) \
                            - (np.sin(de_gp) * cos_lat * cos_lcp_lon))
    dec = np.arcsin(sin_d)
    ra = (ramragp + ra_gp + (2. * np.pi)) % (2. * np.pi)
    return np.degrees(ra), np.degrees(dec)

############################################################

def dec2hms(dec):
    """Convert from decimal degrees to hours-minutes-seconds.

    ADW: This should be replaced by astropy...
    from astropy.coordinates import Angle
    hms = Angle(dec*u.deg).hms
    return (hms.h,hms.m,hms.s)

    """
    DEGREE = 360.
    HOUR = 24.
    MINUTE = 60.
    SECOND = 3600.
    
    dec = float(dec)
    fhour = dec*(HOUR/DEGREE)
    hour = int(fhour)

    fminute = (fhour - hour)*MINUTE
    minute = int(fminute)
    
    second = (fminute - minute)*MINUTE
    return (hour, minute, second)

def dec2dms(dec):
    """Convert from decimal degrees to degrees-minutes-seconds.

    ADW: This should be replaced by astropy
    from astropy.coordinates import Angle
    dms = Angle(dec*u.deg).dms
    return (dms.d,dms.m,dms.s)

    """
    DEGREE = 360.
    HOUR = 24.
    MINUTE = 60.
    SECOND = 3600.

    dec = float(dec)
    sign = np.copysign(1.0,dec)

    fdeg = np.abs(dec)
    deg = int(fdeg)
    
    fminute = (fdeg - deg)*MINUTE
    minute = int(fminute)
    
    second = (fminute - minute)*MINUTE

    # Careful, need float to allow negative zeros
    deg = sign*int(deg)
    return (deg, minute, second)

def hms2dec(hms):
    """Convert longitude from hours,minutes,seconds in string or 3-array
    format to decimal degrees.

    ADW: This really should be replaced by astropy

    """
    DEGREE = 360.
    HOUR = 24.
    MINUTE = 60.
    SECOND = 3600.

    if isstring(hms):
        hour,minute,second = np.array(re.split('[hms]',hms))[:3].astype(float)
    else:
        hour,minute,second = hms.T

    decimal = (hour + minute * 1./MINUTE + second * 1./SECOND)*(DEGREE/HOUR)
    return decimal

def dms2dec(dms):
    """Convert latitude from degrees,minutes,seconds in string or 3-array
    format to decimal degrees.

    ADW: This really should be replaced by astropy
    """
    DEGREE = 360.
    HOUR = 24.
    MINUTE = 60.
    SECOND = 3600.

    # Be careful here, degree needs to be a float so that negative zero
    # can have its signbit set:
    # http://docs.scipy.org/doc/numpy-1.7.0/reference/c-api.coremath.html#NPY_NZERO

    if isstring(dms):
        degree,minute,second = np.array(re.split('[dms]',dms))[:3].astype(float)
    else:
        degree,minute,second = dms.T

    sign = np.copysign(1.0,degree)
    decimal = np.abs(degree) + minute * 1./MINUTE + second * 1./SECOND
    decimal *= sign
    return decimal

def sr2deg(solid_angle):
    """ Convert solid angle from steradians to square deg."""
    return np.degrees(np.degrees(solid_angle))

def deg2sr(solid_angle):
    """ Convert solid angle from square deg to steradians"""
    return np.radians(np.radians(solid_angle))

############################################################

def distanceToDistanceModulus(distance):
    """ Return distance modulus for a given distance (kpc).

    Parameters
    ----------
    distance : distance (kpc)

    Returns
    -------
    mod : distance modulus
    """
    return 5. * (np.log10(np.array(distance)) + 2.)

dist2mod = distanceToDistanceModulus

def distanceModulusToDistance(distance_modulus):
    """ Return distance (kpc) for a given distance modulus.

    Parameters
    ----------
    distance_modulus : distance modulus

    Returns
    -------
    distance : distance (kpc)
    """
    return 10**((0.2 * np.array(distance_modulus)) - 2.)

mod2dist = distanceModulusToDistance

############################################################

def ang2const(lon,lat,coord='gal'):
    import ephem

    scalar = np.isscalar(lon)
    lon = np.array(lon,copy=False,ndmin=1)
    lat = np.array(lat,copy=False,ndmin=1)

    if coord.lower() == 'cel':
        ra,dec = lon,lat
    elif coord.lower() == 'gal':
        ra,dec = gal2cel(lon,lat)
    else:
        msg = "Unrecognized coordinate"
        raise Exception(msg)

    x,y = np.radians([ra,dec])
    const = [ephem.constellation(coord) for coord in zip(x,y)]
    if scalar: return const[0]
    return const

def ang2iau(lon,lat,coord='gal'):
    """
    Convert from coordinates to IAU naming convention.
    Naming has precision of one minute: J{HH}{MM}+{DD}{MM}

    See:
    https://www.iau.org/public/themes/naming/
    http://cdsweb.u-strasbg.fr/Dic/iau-spec.html

    Parameters
    ----------
    lon   : longitude (deg)
    lat   : latitude (deg)
    coord : coordinate system for lon/lat ['gal','cel']

    Returns
    -------
    name  : name consistent with IAU convention
    """
    # Default name formatting
    fmt = "J%(hour)02i%(hmin)02i%(deg)+03.0f%(dmin)02i"

    scalar = np.isscalar(lon)
    lon = np.array(lon,copy=False,ndmin=1)
    lat = np.array(lat,copy=False,ndmin=1)

    if coord.lower() == 'cel':
        ra,dec = lon,lat
    elif coord.lower() == 'gal':
        ra,dec = gal2cel(lon,lat)
    else:
        msg = "Unrecognized coordinate"
        raise Exception(msg)

    iau = []
    for _ra,_dec in zip(ra,dec):
        hms = dec2hms(_ra); dms = dec2dms(_dec)
        params = dict(hour=hms[0],hmin=hms[1],
                      deg=dms[0],dmin=dms[1])
        iau.append(fmt%params)
    if scalar: return iau[0]
    return np.array(iau)


def match(lon1, lat1, lon2, lat2, tol=None, nnearest=1):
    """
    Adapted from Eric Tollerud.
    Finds matches in one catalog to another.
 
    Parameters
    lon1 : array-like
        Longitude of the first catalog (degrees)
    lat1 : array-like
        Latitude of the first catalog (shape of array must match `lon1`)
    lon2 : array-like
        Longitude of the second catalog
    lat2 : array-like
        Latitude of the second catalog (shape of array must match `lon2`)
    tol : float or None, optional
        Proximity (degrees) of a match to count as a match.  If None,
        all nearest neighbors for the first catalog will be returned.
    nnearest : int, optional
        The nth neighbor to find.  E.g., 1 for the nearest nearby, 2 for the
        second nearest neighbor, etc.  Particularly useful if you want to get
        the nearest *non-self* neighbor of a catalog.  To do this, use:
        ``spherematch(lon, lat, lon, lat, nnearest=2)``
 
    Returns
    -------
    idx1 : int array
        Indices into the first catalog of the matches. Will never be
        larger than `lon1`/`lat1`.
    idx2 : int array
        Indices into the second catalog of the matches. Will never be
        larger than `lon2`/`lat2`.
    ds : float array
        Distance (in degrees) between the matches
    """
    from scipy.spatial import cKDTree
 
    lon1 = np.asarray(lon1)
    lat1 = np.asarray(lat1)
    lon2 = np.asarray(lon2)
    lat2 = np.asarray(lat2)
 
    if lon1.shape != lat1.shape:
        raise ValueError('lon1 and lat1 do not match!')
    if lon2.shape != lat2.shape:
        raise ValueError('lon2 and lat2 do not match!')

    rotator = SphericalRotator(0,0)

 
    # This is equivalent, but faster than just doing np.array([x1, y1, z1]).T
    x1, y1, z1 = rotator.cartesian(lon1.ravel(),lat1.ravel())
    coords1 = np.empty((x1.size, 3))
    coords1[:, 0] = x1
    coords1[:, 1] = y1
    coords1[:, 2] = z1
 
    x2, y2, z2 = rotator.cartesian(lon2.ravel(),lat2.ravel())
    coords2 = np.empty((x2.size, 3))
    coords2[:, 0] = x2
    coords2[:, 1] = y2
    coords2[:, 2] = z2
 
    tree = cKDTree(coords2)
    if nnearest == 1:
        idxs2 = tree.query(coords1)[1]
    elif nnearest > 1:
        idxs2 = tree.query(coords1, nnearest)[1][:, -1]
    else:
        raise ValueError('invalid nnearest ' + str(nnearest))
 
    ds = angsep(lon1, lat1, lon2[idxs2], lat2[idxs2])
 
    idxs1 = np.arange(lon1.size)
 
    if tol is not None:
        msk = ds < tol
        idxs1 = idxs1[msk]
        idxs2 = idxs2[msk]
        ds = ds[msk]
 
    return idxs1, idxs2, ds