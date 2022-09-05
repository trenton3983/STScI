import numpy, sys, json, requests, warnings

from astropy.table import Table
from astropy.io import fits
from astropy.wcs import WCS, FITSFixedWarning
from astropy.visualization import PercentileInterval, AsinhStretch

from PIL import Image
from io import BytesIO
import matplotlib.pyplot as plt

from urllib.parse import quote as urlencode
from urllib.request import urlretrieve
import http.client as httplib 


def getcolorim(ra, dec, size=240, zoom=1.0, output_size=None, filters="grizy", format="jpg", **kw):
    
    """Get color image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    zoom = scale factor to reduce image size (default = 1, ignored if output_size is specified).
           E.g., use zoom=0.5 to reduce image size by a factor of 2 in x and y.
           zoom is rounded to the nearest inverse integer value (e.g., zoom=0.3 is the same as
           zoom = 0.3333333333).
           Unlike output_size, zoom can be used for fits format images.
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png")
    Returns the image (a PIL Image)
    """
    
    if format not in ("jpg","png"):
        raise ValueError("format must be jpg or png")
    url = geturl(ra, dec, size=size, zoom=zoom, output_size=output_size,
                 filters=filters, format=format, color=True, **kw)
    r = requests.get(url)
    im = Image.open(BytesIO(r.content))
    return im


def getgrayim(ra, dec, size=240, zoom=1.0, output_size=None, filter="g", format="jpg", **kw):
    
    """Get grayscale image at a sky position
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    zoom = scale factor to reduce image size (default = 1, ignored if output_size is specified).
           E.g., use zoom=0.5 to reduce image size by a factor of 2 in x and y.
           zoom is rounded to the nearest inverse integer value (e.g., zoom=0.3 is the same as
           zoom = 0.3333333333).
           Unlike output_size, zoom can be used for fits format images.
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filter = string with filter to extract (one of grizy)
    format = data format (options are "jpg", "png" or "fits")
    Returns the image (either a FITS hdu or a PIL Image)
    """
    
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be jpg, png or fits")
    if filter not in list("grizy"):
        raise ValueError("filter must be one of grizy")
    url = geturl(ra, dec, size=size, zoom=zoom, output_size=output_size,
                 filters=filter, format=format, **kw)
    r = requests.get(url[0])
    if format == "fits":
        im = fits.open(BytesIO(r.content))
    else:
        im = Image.open(BytesIO(r.content))
    return im


def getwcs(image, extension=0, verbose=False):
    """Get WCS for either a FITS image or a PIL JPEG Image
    
    Note that only JPEG images from MAST fitscut cutout services (such as 
    ps1images.stsci.edu, archive.stsci.edu and hla.stsci.edu) will have the necessary
    WCS information.
    """
    # try FITS first
    try:
        # suppress minor warnings
        with warnings.catch_warnings():
            warnings.simplefilter('ignore', FITSFixedWarning)
            return WCS(image[extension].header)
    except (TypeError, AttributeError):
        pass
    if image.format != 'JPEG':
        raise ValueError(f"cannot get WCS from PIL image of type '{image.format}'")
    wstring = image.app.get('COM')
    if not wstring:
        raise ValueError("no WCS data found in JPEG comment")
    if verbose:
        print("wcs string")
        print(wstring.decode())
    w = WCS(fits.Header.fromstring(wstring,sep='\n'))
    return w


def getimages(ra,dec,filters="grizy",verbose=False):
    
    """Query ps1filenames.py service to get a list of images
    
    ra, dec = position in degrees
    size = image size in pixels (0.25 arcsec/pixel)
    filters = string with filters to include
    Returns a table with the results
    """
    
    service = "https://ps1images.stsci.edu/cgi-bin/ps1filenames.py"
    url = f"{service}?ra={ra}&dec={dec}&filters={filters}"
    table = Table.read(url, format='ascii')
    if verbose:
        print(f"Found {len(table)} images for ra {ra} dec {dec} filters '{filters}'")
    return table


def geturl(ra, dec, size=240, zoom=1.0, output_size=None, filters="grizy", format="jpg", color=False, autoscale=99.75, asinh=True, verbose=False):
    
    """Get URL for images in the table
    
    ra, dec = position in degrees
    size = extracted image size in pixels (0.25 arcsec/pixel)
    zoom = scale factor to reduce image size (default = 1, ignored if output_size is specified).
           E.g., use zoom=0.5 to reduce image size by a factor of 2 in x and y.
           zoom is rounded to the nearest inverse integer value (e.g., zoom=0.3 is the same as
           zoom = 0.3333333333).
           Unlike output_size, zoom can be used for fits format images.
    output_size = output (display) image size in pixels (default = size).
                  output_size has no effect for fits format images.
    filters = string with filters to include
    format = data format (options are "jpg", "png" or "fits")
    color = if True, creates a color image (only for jpg or png format).
            Default is return a list of URLs for single-filter grayscale images.
    autoscale and asinh set image contrast for jpg and png format (ignored for fits)
    Returns a string with the URL for a color image or a list of URLs for grayscale images
    """
    
    if color and format == "fits":
        raise ValueError("color images are available only for jpg or png formats")
    if format not in ("jpg","png","fits"):
        raise ValueError("format must be one of jpg, png, fits")
    table = getimages(ra,dec,filters=filters,verbose=verbose)
    url = "https://ps1images.stsci.edu/cgi-bin/fitscut.cgi?"
    params = [ f"ra={ra}", f"dec={dec}", f"format={format}" ]
    if output_size and format != "fits":
        params.append(f"output_size={output_size}")
    elif zoom != 1:
        params.append(f"zoom={zoom}")
        # round size up to a multiple of 1/zoom (not strictly necessary, but may
        # improve quality of edge pixels)
        zround = round(1./zoom)
        if zround > 1:
            size += (zround-size) % zround
    params.append(f"size={size}")
    if format != "fits":
        params.append(f"autoscale={autoscale}")
        if asinh:
            params.append(f"asinh=1")
    # sort filters from red to blue
    flist = ["yzirg".find(x) for x in table['filter']]
    table = table[numpy.argsort(flist)]
    if color:
        if len(table) > 3:
            # pick 3 filters
            table = table[[0,len(table)//2,len(table)-1]]
        for i, param in enumerate(["red","green","blue"]):
            params.append("{}={}".format(param,table['filename'][i]))
        url = url + "&".join(params)
    else:
        params.append("red=")
        urlbase = url + "&".join(params)
        url = []
        for filename in table['filename']:
            url.append(urlbase+filename)
    if verbose:
        print("url", url)
    return url
    
    
def getgaia(ra, dec, radius, epoch=None, version="edr3", verbose=False, url='https://gsss.stsci.edu/webservices/VO/CatalogSearch.aspx'):

    """Get Gaia sources at the given position

    Returns table of Gaia sources.  If epoch is specified, new fields at
    that epoch (nra, ndec) are computed using the proper motions, along
    with errors (nra_error, ndec_error, nra_dec_corr) computed using the
    Gaia errors and covariances.  If the epoch is omitted, no PM
    corrections are made.

    :param ra:         Search RA, units in degrees.
    :param dec:        Search Dec, units in degrees.
    :param radius:     Search radius, units in degrees.
    :param epoch:      Epoch to apply PMs (default is reference epoch)
    :param version:    Version of the Gaia catalog. 'dr2' or 'edr3'.
    :param verbose:    If true, prints some info
    :param url:        URL for Gaia server
    :type ra:          float
    :type dec:         float
    :type radius:      float
    :type version:     str
    :type verbose:     bool
    :type url:         str
    :returns:          (tab) astropy.Table: Gaia catalog table with (ra,dec) coordinates, magnitudes,
                       errors, etc., shifted to epoch (if specified)
    """
    vlist = ['dr2', 'edr3']
    if version not in vlist:
        raise ValueError("version '{}' must be one of {}".format(version,', '.join(vlist)))
    catname = 'gaia'+version
    r = requests.get(url, params=dict(ra=ra, dec=dec, sr=radius, version=version,
            format='csv', cat=catname))
    try:
        gcat = Table.read(r.text,format='ascii.csv',comment='#')
        # change column names to lowercase
        for col in gcat.colnames:
            lcol = col.lower()
            if lcol != col:
                gcat.rename_column(col,lcol)
                if verbose:
                    print("renamed {} to {}".format(col,lcol))
    except IndexError:
        # no Gaia sources found
        gcat = []
    except Exception as e:
        print(r.text)
        raise

    if verbose:
        print(f"got {len(gcat)} Gaia sources in cone")
    if len(gcat) == 0:
        return gcat
    if epoch:
        applypm(gcat, epoch)
    return gcat


def applypm(tab, epoch):
    """
    Apply the proper motions between the Gaia catalog epoch and the specified
    observing epoch.  Modify tab by adding new columns nra, ndec, nra_error,
    ndec_error, nra_dec_corr and nepoch with the positions shifted to the
    target epoch using the proper motions (but not the parallax at the
    moment).  The PM deltas (in mas) ndra,nddec are also added.  The errors
    are modified for the epoch using the catalog errors and covariances.

    The units for the new quantities are the same as for the original
    catalog values (deg for nra and ndec, mas for the errors, ndra and ndec,
    and years for nepoch).

    The positions will have null values for objects without known proper
    motions.

    tab is modified in place by adding the new columns.

    :param tab:        Gaia catalog table
    :param epoch:      Desired epoch for positions (decimal years)
    :type tab:         (tab) astropy.Table: Gaia catalog table with (ra,dec) coordinates, magnitudes etc.
    :type epoch:       float
    :returns:          None
    """
    dt = epoch - tab['ref_epoch']
    dra = tab['pmra']*dt
    ddec = tab['pmdec']*dt
    tab['ndra'] = dra
    tab['nddec'] = ddec
    tab['nra'] = tab['ra'] + dra*(1.e-3/(3600.0*np.cos((np.pi/180)*tab['dec'])))
    tab['ndec'] = tab['dec'] + ddec*(1.e-3/3600.0)
    ra_var = (tab['ra_error']**2 +
              (dt*tab['pmra_error'])**2 +
              2*dt*tab['ra_pmra_corr']*tab['ra_error']*tab['pmra_error'])
    dec_var = (tab['dec_error']**2 +
               (dt*tab['pmdec_error'])**2 +
               2*dt*tab['dec_pmdec_corr']*tab['dec_error']*tab['pmdec_error'])
    ra_dec_covar = (tab['ra_dec_corr']*tab['ra_error']*tab['dec_error'] +
                    dt*tab['ra_pmdec_corr']*tab['ra_error']*tab['pmdec_error'] +
                    dt*tab['dec_pmra_corr']*tab['dec_error']*tab['pmra_error'] +
                    dt**2*tab['pmra_pmdec_corr']*tab['pmra_error']*tab['pmdec_error'])
    tab['nra_error'] = np.sqrt(ra_var)
    tab['ndec_error'] = np.sqrt(dec_var)
    tab['nra_dec_corr'] = ra_dec_covar/np.sqrt(ra_var*dec_var)
    tab['nepoch'] = epoch
    for c in ['nra_error', 'ndec_error', 'ndra', 'nddec',
            'ra_error','dec_error','pmra','pmdec','pmra_error','pmdec_error']:
        tab[c].format = ".2f"
    tab['nra_dec_corr'].format = ".4f"
    tab['nepoch'].format = ".2f"
    tab['nra'].format = ".12f"
    tab['ndec'].format = ".13f"