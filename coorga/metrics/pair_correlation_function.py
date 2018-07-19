#!/usr/bin/env python


import numpy as np
import scipy.integrate, scipy.interpolate, scipy.spatial


import tropy.analysis_tools.grid_and_interpolation as gi



######################################################################
######################################################################

##LLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLL
## automatically generated by 10-PairCorrelationFunction.ipynb
##TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT

def pairCorrelationFunction_2D(pos, egrid, numberDensity, rMax, dr, 
                               rbins = None,
                               equal_area = True, 
                               weight = None,
                               constant_analytic = False,
                               mask_edge_particle = True,
                               use_radial_interpolation = False,
                               finite_size_correction = False):

    """Compute the two-dimensional pair correlation function, also known
    as the radial distribution function, for a set of circular particles
    contained in a square region of a plane.  

    This simple function finds reference particles such that a circle of radius rMax drawn around the
    particle will fit entirely within the square, eliminating the need to
    compensate for edge effects.  If no such particles exist, an error is
    returned. Try a smaller rMax...or write some code to handle edge effects! ;)


    Parameters
    ----------
        
    pos : numpy array
        position vector such that x, y = pos.T and pos.shape = (Nparticles, 2)
        
    egrid : numpy array
        edge-based grid at which number density field is given
    
    numberDensity : numpy array
        number density field, as number of particles per area
        
    rMax : float value
        maximum radius
 
    dr : float value
        increment for increasing radius of annulus
        
    rbins : numpy array, optional, default = None
        array which contains range bin values

    equal_area : bool, optional, default = True
        switch that determines if equi-distance radius rings or rings of equal area
        are chosen, optional
        
    weight : float, optional, default = None
        if weight (number between 0 and 1) is set, than radius rings are calculated 
        from a weighted version between equi-distant (weight = 0) and equal area (weight = 1)


    constant_analytic : bool, optional default = False
        use analytic form of expected reference number assuming a constant number density
 
    mask_edge_particles : bool, optional default = True
        If particle close to edge are ignored
 
    use_radial_interpolation : bool, optional,  default = False
        if number field is represented by an interpolating function and intergation is performed
        on that field

    finite_size_correction : bool, optional, default = False
        if a finite size correction N / (N-1) is applied


    Returns 
    -------
    
    g(r) : numpy array 
        correlation function g(r) for each particle
        
    g_ave(r) : numpy array 
        average correlation function
        
    radii : numpy array 
        radii of the   annuli used to compute g(r)

     interior_indices : numpy array
        indices of reference particles
        
    Note
    ------
    Possibly outdated! Try the faster variant pcf that uses precalculated reference number fields!
    
    """
    # Number of particles in ring/area of ring/number of reference particles/number density
    # area of ring = pi*(r_outer**2 - r_inner**2)

    # Find particles which are close enough to the box center that a circle of radius
    # rMax will not cross any edge of the box

    # get particle coordiantes
    # =========================
    x, y = pos.T 
    
    
    # get grids
    # ==========
    xg, yg = egrid   # edge based grid
    xc = gi.lmean(gi.lmean(xg, axis = 0), axis = 1)   # center point based grid
    yc = gi.lmean(gi.lmean(yg, axis = 0), axis = 1)   # center point based grid
    
    
    # check boundaries
    # ================

    # overwrite rMax if rbins is given!!!
    if rbins != None:
        rMax = rbins.max()
    
    xmin, xmax = xg.min(), xg.max()
    ymin, ymax = yg.min(), yg.max()
    
    if mask_edge_particle:
        xmask = (x > rMax + xmin) & (x < xmax - rMax)
        ymask = (y > rMax + ymin) & (y < ymax - rMax)
    else:
        xmask = (x > xmin) & (x < xmax)
        ymask = (y > ymin) & (y < ymax)
    
    interior_indices, = np.where( xmask & ymask )
    num_interior_particles = len(interior_indices)

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a circle of radius rMax\
                will lie entirely within a square of side length S.  Decrease rMax\
                or increase the size of the square.")

        
    # count particles per distance ring
    # ==================================
    if rbins == None:
        edges = radius_ring_edges(rMax, dr, equal_area = equal_area, weight = weight)
    else:
        edges = rbins
        
    num_increments = len(edges) - 1
    g = np.zeros([num_interior_particles, num_increments])

    # finite size correction (N - 1) / N
    # =================================
    ncells = 1. * pos.shape[0]
    if finite_size_correction:
        corr_factor = (ncells - 1.) / ncells
    else:
        corr_factor = 1.
    
    # Compute pairwise correlation for each interior particle
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = np.sqrt((x[index] - x)**2 + (y[index] - y)**2)
        d[index] = 2 * rMax

        (npart, bins) = np.histogram(d, bins = edges, normed=False)
        
        dgrid = np.sqrt((x[index] - xc)**2 + (y[index] - yc)**2)
        
        # here, we calculate reference number of particles
        for i in range(num_increments):        
            r1 = edges[i]
            r2 = edges[i + 1]
            
            if constant_analytic:
                n0 = np.pi * numberDensity.mean() * (r2**2 - r1**2)
            else:
                if use_radial_interpolation:
                    p0 = (x[index], y[index])
                    n0 = radial_ring_integration_2d_field(egrid, numberDensity, r1, r2, p0 = p0)
                else:
                    dmask = (dgrid >= edges[i]) & (dgrid <= edges[i + 1])
                    n0 = integrate_2dfield(egrid, numberDensity, dmask)
            
            g[p, i] = np.ma.divide( npart[i] ,  n0 )
            

    # Average g(r) for all interior particles and compute radii
    #g_average = zeros(num_increments)
    #for i in range(num_increments):
    #    radii[i] = (edges[i] + edges[i+1]) / 2.
    #    rOuter = edges[i + 1]
    #    rInner = edges[i]
    #    g_average[i] = mean(g[:, i]) / (pi * (rOuter**2 - rInner**2))

    g = corr_factor * np.ma.masked_invalid( g )
    
    g_average = g.mean( axis = 0)
    radii = gi.lmean(edges)

    return (g, g_average, radii, interior_indices)



######################################################################
######################################################################



def pcf(pos, rbins, domain, refgrid, ref_number, tree_as_grid_input = False):
    
    """
    Compute the two-dimensional pair correlation function, also known
    as the radial distribution function, for a set of circular particles
    contained in a square region of a plane.  

    Parameters
    ----------
   
    pos : numpy array
        position vector such that x, y = pos.T and pos.shape = (Nparticles, 2)
        
    rbins : numpy array
        ring edge vector
    
    domain : tuple or list in the form ((xmin, xmax), (ymin, ymax))
        containing min & max coordiantes for masking the edge
        
    refgrid :  tuble of two 2dim numpy arrays OR scipy tree class instance
        grid on which reference number of cells is given
        
    ref_number : numpy array
        reference number of cells per ring

    tree_as_grid_input : bool, optional, default = False
        if not reference grid, but scipy tree class instance (for nearest neighbor int) is supplied


    Returns 
    -------
    
    g(r) : numpy array 
        correlation function g(r) for each particle
        
    g_ave(r) : numpy array 
        average correlation function
        
    radii : numpy array 
        radii of the   annuli used to compute g(r)

     interior_indices : numpy array
        indices of reference particles
    """

    # Number of particles in ring/area of ring/number of reference particles/number density
    # area of ring = pi*(r_outer**2 - r_inner**2)

    # Find particles which are close enough to the box center that a circle of radius
    # rMax will not cross any edge of the box

    # get particle coordiantes
    # =========================
    x, y = pos.T 
    
    
    
    
    # check boundaries
    # ================
    xmin, xmax = domain[0]
    ymin, ymax = domain[1]
        
    rMax = rbins.max()
    
    xmask = (x > rMax + xmin) & (x < xmax - rMax)
    ymask = (y > rMax + ymin) & (y < ymax - rMax)
    
    interior_indices, = np.where( xmask & ymask )
    num_interior_particles = len(interior_indices)
    
    pos_inner = pos[interior_indices]

    if num_interior_particles < 1:
        raise  RuntimeError ("No particles found for which a circle of radius rMax\
                will lie entirely within a square of side length S.  Decrease rMax\
                or increase the size of the square.")

        
    # set distance rings
    # ==================
    edges = rbins
    num_increments = len(edges) - 1

    g = np.zeros([num_interior_particles, num_increments])
    npart = np.zeros([num_interior_particles, num_increments])

    
    # Compute number of observed particles
    # ====================================
    for p in range(num_interior_particles):
        index = interior_indices[p]
        d = np.sqrt((x[index] - x)**2 + (y[index] - y)**2)
        d[index] = 2 * rMax

        (npart[p,:], bins) = np.histogram(d, bins = edges, normed = False)
        

    # Get number of reference particles
    # =================================
    if tree_as_grid_input:
        grid_tree = refgrid
    else:
        xout, yout = refgrid
        pout = np.array([xout.flatten(), yout.flatten()]).T
        grid_tree = scipy.spatial.KDTree(pout)
    
    dist, index = grid_tree.query(pos_inner)
    
    nref = ref_number.reshape(-1, num_increments)[index]

    
    g = np.ma.divide( npart, nref )
    g = np.ma.masked_invalid( g )
    
    g_average = g.mean( axis = 0)
    radii = gi.lmean(edges)

    return (g, g_average, radii, interior_indices)


######################################################################
######################################################################


def init_reference_numbers(egrid, numberDensity, rbins, nsub = 4, radial_integration = True):
    
    '''
    Calculates the number of expected cells / particle as function of distance.
    

    Parameters
    -----------
    egrid : tuble of list of two 2dim numpy arrays
        edge-based grid at which number density field is given
        
    numberDensity : numpy array
        number density field, as number of particles per area
        
    rbins : numpy array
        ring edge array

    nsub : int values, optional, default = 4
        values that determines subsampling of the analyzed fields

    radial_integration : bool, optional, default = True
        switch if radial integration of an interpolated representation of the 
        number field is used.

        
    Returns
    --------
    xout : numpy array, 2dim with shape = (nrows, ncols)
       output x-grid
 
    yout : numpy array, 2dim with shape = (nrows, ncols)
       output y-grid
    
    n0 : numpy array, 3dim with shape = (nrows, ncols, nrbins)
       expected number of cells on grid for each ring interval
    '''

    
    # get grids
    # ==========
    xg, yg = egrid   # edge based grid
    xc = gi.lmean(gi.lmean(xg, axis = 0), axis = 1)   # center point based grid
    yc = gi.lmean(gi.lmean(yg, axis = 0), axis = 1)   # center point based grid
    
    dx = gi.lmean( np.diff(xg, axis = 0), axis = 1)
    dy = gi.lmean( np.diff(yg, axis = 1), axis = 0)
       
    n_dx_dy = numberDensity * dx * dy
        
    
    # set output grid
    # =================
    xout = xc[::nsub, ::nsub]
    yout = yc[::nsub, ::nsub]
    
    pc = np.array([xc.flatten(), yc.flatten()]).T
    grid_tree = scipy.spatial.KDTree(pc)

    nrow_out, ncol_out = xout.shape
    
    
    # count particles per distance ring
    # ==================================
    edges = rbins
    num_increments = len(edges) - 1
    
    n0 = np.zeros((nrow_out, ncol_out, num_increments))
    ncum = np.zeros((nrow_out, ncol_out, num_increments))
    
    
    # calculated expected number of reference cells
    # ==============================================
    for irow in range(nrow_out):
        for icol in range(ncol_out):
            
            p0 = (xout[irow, icol], yout[irow, icol])
            
            if radial_integration:
                n0[irow, icol, :] = radial_integration_2d_field(egrid, numberDensity, rbins, p0)
            else:

                for ibin, r in enumerate(edges[1:]):
                    m = grid_tree.query_ball_point(p0, r)
                    ncum[irow, icol, ibin] = n_dx_dy.flatten()[m].sum()

                n0[irow, icol, 0] = ncum[irow, icol, 0]
                n0[irow, icol, 1:] = np.diff( ncum[irow, icol, ::-1] )
    return xout, yout, n0


#################################################################
#################################################################


def integrate_2dfield(egrid, f, mask):
    
    '''
    Integrates a 2d field over a certain region given a mask and grid.
    
    
    Parameters
    ----------
    egrid : tuble of list of numpy arrays
        edge-based grid such that:  xeg, yeg = egrid

    f : numpy array, 2dim
        input field that is integrated

    mask : numpy array, 2dim, bool
        mask of a certain region that is used in integration (same shape as f)
    

    Returns
    -------
    fint : float value
        integral of f for the region where mask is True
    '''
    
    xeg, yeg = egrid
    
    dx = gi.lmean( np.diff(xeg, axis = 0), axis = 1)
    dy = gi.lmean( np.diff(yeg, axis = 1), axis = 0)
       

    fint = (f[mask] * dx[mask] * dy[mask]).sum()
    
    return fint


##########################################################################
##########################################################################

def transform2cylinder_coords(xm, ym, f, p0 = (0,0)):
    
    '''
    Interpolates and transforms a function given in Cartesian coordinates (x,y)
    into cylindrical coordinates (r, phi).
    
    Parameters
    ------------
    xm : numpy array, 2dim
        center-based x-grid
    
    ym : numpy array, 2dim 
        center-based y-grid

    f : numpy array, 2dim
        gridded field

    p0 : tuple of 2 float values, optional, default = (0,0)
        base point, i.e. x0, y0 = p0
    

    Returns
    --------
    ftrans: function of two arguments (r, phi)
        interpolation function in (r, phi)-coordinates
    '''
    
    # get base point
    x0, y0 = p0
    
    # get relative coordinates
    dx = xm - x0
    dy = ym - y0
    
    # interpolation
    fint = scipy.interpolate.RectBivariateSpline(dx, dy, f)
    
    
    # coordinate transformation
    x = lambda r, phi: r * np.cos( phi )
    y = lambda r, phi: r * np.sin( phi )
    
    ftrans = lambda r, phi: fint.ev( x(r, phi), y(r, phi) )
    
    return ftrans

##########################################################################
##########################################################################

def radial_ring_integration_2d_field(egrid, f, r1, r2, p0):
    
    '''
    Integrates a 2d field over a certain region given a base point 
    and the ring edges.
    
    
    Parameters
    ----------
    egrid : tuple or list of two numpy arrays
        edge-based grid such that:  xeg, yeg = egrid

    f : numpy array, 2dim
        2d field

    r1 : float value
        inner ring radius

    r2 : float value
        outer ring radius

    p0 : tuple of 2 float values
        base point, i.e. x0, y0 = p0


    Returns
    --------
    fint : float value
        radial integral of f
    '''
    
    
    # get the equi-distant grids
    xeg, yeg = egrid
    
    xe = xeg.mean(axis = 1)
    ye = yeg.mean(axis = 0)
    
    xm = gi.lmean( xe )    
    ym = gi.lmean( ye )
    
   
    # transformed function: now in cylinder coordinates
    ftrans = transform2cylinder_coords(xm, ym, f, p0 = p0)
    kernel = lambda r, phi: r * ftrans(r, phi)
    
    I = scipy.integrate.nquad(kernel, [[r1, r2], [0, 2*np.pi]])
    return I[0]
    

##########################################################################
##########################################################################

def radial_integration_2d_field(egrid, f, rbins, p0):
    
    '''
    Integrates a 2d field over a certain region given a base point 
    and the ring edges.
    
    
    Parameters
    ----------
    egrid :  tuple or list of two numpy arrays
        edge-based grid such that:  xeg, yeg = egrid

    f : numpy array, 2dim
        2d field

    rbins : numpy array
       array of range bins

    p0 : tuple of 2 float values
        base point, i.e. x0, y0 = p0

    
    Returns
    --------
    fint : float
        radial integral of f
    '''
    
    
    # get the equi-distant grids
    xeg, yeg = egrid
    
    xe = xeg.mean(axis = 1)
    ye = yeg.mean(axis = 0)
    
    xm = gi.lmean( xe )    
    ym = gi.lmean( ye )
    
   
    # transformed function: now in cylinder coordinates
    ftrans = transform2cylinder_coords(xm, ym, f, p0 = p0)
    kernel = lambda r, phi: r * ftrans(r, phi)
    
    nbins = rbins.shape[0]
    I = np.zeros( nbins - 1 )
    for i in range(nbins - 1):
        r1 = rbins[i]
        r2 = rbins[i + 1]
        I[i] = scipy.integrate.nquad(kernel, [[r1, r2], [0, 2*np.pi]])[0]
    return I
    
#################################################################
#################################################################



def radius_ring_edges(rMax, dr, equal_area = True, weight = None):

    '''
    Returns inner and outer ring edges for circular ring elements.
    

    Parameters
    -----------
    rMax : float value
        maximum ring radius

    dr : float value
        ring element size

    equal_area :  bool, optional, default = True
        if True, ring elements are calculated to have equal area

    weight: float, optional, default = None
        if number between 0 and 1, a weighted version between 
        equidistant and equal-area rings is chosen
            

    Returns
    --------
    redges : numpy array
        ring edges
    
    '''
    # weight variant
    if weight != None:
        r_A_equal = radius_ring_edges(rMax, dr, equal_area = True, weight = None)
        r_equi_dist = radius_ring_edges(rMax, dr, equal_area = False, weight = None)
        
        redges = weight * r_equi_dist + (1 - weight) * r_A_equal
        
    
    if equal_area and weight == None:
        
        # largest area
        A_max = np.pi * rMax**2

        # individual ring area
        nrings = np.int(rMax / dr)
        A_ring = A_max /  nrings

        # diamter of ring-equivalent circle
        rq = A_ring / (np.pi)

        # iteratively calculate outer ring edges / first inner ring is at zero -> inner circle
        redges = np.zeros(nrings + 1)

        for i in range(nrings):
            redges[i + 1] = np.sqrt(rq + redges[i]**2)
            
    elif not equal_area and weight == None: # make equal distance
        redges = np.arange(0., rMax + 0.1 * dr, dr)
        
        
    
    return redges

######################################################################
######################################################################
