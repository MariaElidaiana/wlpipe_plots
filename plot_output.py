#!/usr/env python
"""
.. module:: plot_output
:synopsis: Script to plot the data vector results from WLpipe and compare with DES Y1.
.. moduleauthor:: Maria Elidaiana <mariaeli@brandeis.edu>
"""

import sys
from astropy.table import Table as tbl
from matplotlib import pyplot as plt
import numpy as np

def get_fits(path, filename, hdu):
    """
    Reads the data vector as astropy.table
        
    Parameters
    ----------
    path: str
         Path for the data vector
    filename: str
         Data vector filename
    hdu: int
         hdu number
    
    Returns
    -------
    
    Notes
    -----
    Can be used both to :math:`\\xi_{+}` and :math:`\\xi_{-}`.
    """
    fitstable = tbl.read(path + filename, hdu=hdu)
    fitsbin1 = np.unique(fitstable['BIN1'])
    fitsbin2 = np.unique(fitstable['BIN2']) # The binnig might differ, e.g., for w(theta)
    return fitstable, fitsbin1, fitsbin2

def get_xi_values(dvect, xi_bin, idxs):
    """
    Get the :math:`\\xi` values for each tomographic bin.
    
    Parameters
    ----------
    dvect: data vector
        Data vector
    xi_bin: np.array
        Array with the number of the tomographic bins
    idxs: list of tuples
        List with the tuples (bin1,bin2), where bin1 and bin2 are the tomographic bins range to be used
        
    Returns
    -------
    ang: astropy.table
        Angle in arcmin
    value: astropy.table
        Value of :math:`\\xi`
    rbin: tuple
        Tuple with the tomographic bins range
        
    Notes
    -----
    Can be used both to :math:`\\xi_{+}` and :math:`\\xi_{-}`.
    """
    idx1, idx2 = idxs
    msk   = (dvect['BIN1']==xi_bin[idx1])&(dvect['BIN2']==xi_bin[idx2])
    ang   = dvect['ANG'][msk]
    value = dvect['VALUE'][msk]
    rbin = (xi_bin[idx1], xi_bin[idx2])
    return ang, value, rbin

def get_gammat_values(dvect, idxs):
    """
    Get the :math:`\\gamma_t` values for each tomographic bin.
        
    Parameters
    ----------
    dvect: data vector
        Data vector
    idxs: tuple
        Tuple with the tomographic bins range
        
    Returns
    -------
    ang: astropy.table
        Angles in arcmin
    theta_gammat: astropy.table
        Values of :math:`\\theta\\gamma_t` in units of :math:`10^{-2}` arcmin
    rbin: tuple
        Tuple with the tomographic bins range
    """
    idx1, idx2 = idxs
    msk   = (dvect['BIN1']==idx1)&(dvect['BIN2']==idx2)
    ang   = dvect['ANG'][msk]
    value = dvect['VALUE'][msk]
    theta_gammat = value*ang*1e2
    rbin = (idx1, idx2)
    return ang, theta_gammat, rbin

def get_w_values(dvect, idxs):
    """
    Get the w values for each tomographic bin.
        
    Parameters
    ----------
    dvect:
        Data vector
    idxs: tuple
        Tuple with the tomographic bins range
    
    Returns
    -------
    ang: astropy.table
        Angle in arcmin
    theta_w: astropy.table
        The value of :math:`\\theta w` in units of arcmin
    """
    idx1, idx2 = idxs
    msk   = (dvect['BIN1']==idx1)&(dvect['BIN2']==idx2)
    ang   = dvect['ANG'][msk]
    value = dvect['VALUE'][msk]
    theta_w = value*ang
    rbin = (idx1, idx2)
    return ang, theta_w, rbin

def get_nz_values(dvect):
    """
    Get the n(z) distribution for each tomographic bin of the
    lens and source samples.
    
    Parameters
    ----------
    dvect: data vector
        Data vector
        
    Returns
    -------
    zlow: astropy.table
        Low redshift n(z) distribution
    zmid: astropy.table
        Middle redshift n(z) distribution
    zhigh: astropy.table
        High redshift n(z) distribution
    bin1: astropy.table
        The n(z) of the 1st tomopraphic bin
    bin2: astropy.table
        The n(z) of the 2nd tomopraphic bin
    bin3: astropy.table
        The n(z) of the 3rd tomopraphic bin
    bin4: astropy.table
        The n(z) of the 4th tomopraphic bin
    bin5: astropy.table
        The n(z) of the 5th tomopraphic bin
    
    Notes
    -----
    Only the lens sample has the 5th tomographic bin. For the source
    sample, the function return :code:`None` for the 5th bin.
    """
    zlow  = dvect['Z_LOW']
    zmid  = dvect['Z_MID']
    zhigh = dvect['Z_HIGH']
    bin1  = dvect['BIN1']
    bin2  = dvect['BIN2']
    bin3  = dvect['BIN3']
    bin4  = dvect['BIN4']
    try:
         bin5  = dvect['BIN5']
    except KeyError:
         bin5 = None
    return zlow, zmid, zhigh, bin1, bin2, bin3, bin4, bin5


def plot_xi(dvect1, dvect2, xi_bin, nrow, nidx, bins, ylabel, legend, outname):
    """
    Plot the shear-shear :math:`xi` comparison between the 2pt_pipeline and
    WLpipe on DES Y1 data.
        
    Parameters
    ----------
    dvect1: data vector
        Data vector produced by :code:`2pt_pipeline`
    dvect2: data vector
        Data vector produced by :code:`WLpipe`
    xi_bin: np.array
        Array with the values of the tomographic bins
    nrow: int
        Number of rows for the pyplot.subplots
    nidx: list
        List of indexes for each subplot
    bins: list of tuples
        List of tuples (bin1,bin2) with the tomographic bin combinations
    ylabel: str
        The y-axis label
    legend: str
        Legend for subplots
    outname: str
        Output name for the plot
    
    Notes
    -----
    The y-axis label is to set :math:`\\xi_{+}` and :math:`\\xi_{-}`
    components. The output is a :code:`.png` file with the shear-
    shear comparison between :code:`2pt_pipeline` and :code:`WLpipe`.
    """
    
    plt.figure(1, figsize=(28,22))
    for i in range(len(bins)):
        angle1, value1, rbin = get_xi_values(dvect1, xi_bin, bins[i])
        angle2, value2, rbin = get_xi_values(dvect2, xi_bin, bins[i])
        
        plt.subplot(nrow, nrow, nidx[i])     # two row, two columns, position 1
        if nidx[i]==11:
            plt.axis('off')
            plt.xlim(1e0, 260)
            plt.ylim(1e-8, 1e-3)
            plt.text(1e2, 1e-4, legend, fontsize=25)
            plt.text(1e2, 1e-4, 'blue line=Y1 reference', color='C0', fontsize=25)
        else:
            plt.plot(angle1, value1, lw=2.5)    # y1
            plt.plot(angle2, value2, 'ko', mew=2, ms=12, mfc='None')
            plt.xlim(1e0, 260)
            plt.ylim(1e-8, 1e-3)
            plt.yscale('log', nonposy='mask')
            plt.xscale('log', nonposx='mask')
            plt.text(1e2, 1e-4,'%d, %d'%(rbin[0],rbin[1]), fontsize=25)
            plt.xlabel(r'$\theta$' + ' (arcmin)', fontsize=15)
            plt.ylabel(ylabel, fontsize=15) #(r'$\xi_{+}(\theta)$', fontsize=15)
            plt.tick_params(which='major', length=6, width=1, labelsize=15)
            plt.tick_params(which='minor', length=4, width=1)
    plt.savefig(outname, bbox_inches='tight', dpi=300)
    plt.close()


def plot_gammat(dvect1, dvect2, nrow, ncol, nidx, bins, ylabel, outname):
    """
    Plot the comparison of :math:`\\theta\\gamma_t` (in units of :math:`10^{-2}`
    arcmin) between the :code:`2pt_pipeline` and :code:`WLpipe` on DES Y1 data.
    
    Parameters
    ----------
    dvect1: data vector
        Data vector produced by :code:`2pt_pipeline`
    dvect2: data vector
        Data vector produced by :code:`WLpipe`
    nrow: int
        Number of rows for the pyplot.subplots
    nidx: list
        List of indexes for each subplot
    bins: list of tuples
        List of tuples (bin1,bin2) with the tomographic bin combinations
    ylabel: str
        The y-axis label
    outname: str
        Output name for the plot
    
    Notes
    -----
    The y-axis label is to set :math:`\\xi_{+}` and :math:`\\xi_{-}`
    components. The output is a :code:`.png` file with the shear-
    shear comparison between :code:`2pt_pipeline` and :code:`WLpipe`.
    """
    plt.figure(1, figsize=(38,22))
    for i in range(len(bins)):
        angle1, value1, rbin = get_gammat_values(dvect1, bins[i])
        angle2, value2, rbin = get_gammat_values(dvect2, bins[i])
        
        plt.subplot(nrow, ncol, nidx[i])     # two row, two columns, position 1
        plt.plot(angle1, value1, lw=2.5)    # y1
        plt.plot(angle2, value2, 'ko', mew=2, ms=12, mfc='None')
        plt.axhline(y=0, color='k', linestyle='--')
        plt.xlim(1e0, 450)
        plt.ylim(-0.5, 2.0)
        #plt.yscale('log', nonposy='mask')
        plt.xscale('log', nonposx='mask')
        plt.text(2.0, 1.5,'%d, %d'%(rbin[0],rbin[1]), fontsize=25)
        plt.xlabel(r'$\theta$' + ' (arcmin)', fontsize=15)
        plt.ylabel(ylabel, fontsize=15) #(r'$\xi_{+}(\theta)$', fontsize=15)
        plt.tick_params(which='major', length=6, width=1, labelsize=15)
        plt.tick_params(which='minor', length=4, width=1)
    plt.savefig(outname, bbox_inches='tight', dpi=300)
    plt.close()


def plot_w(dvect1, dvect2, nrow, ncol, nidx, bins, ylabel, legend, outname):
    """
    Plot the comparison of :math:`\\theta w` (in units of arcmin) between the
    :code:`2pt_pipeline` and :code:`WLpipe` on DES Y1 data.
        
    Parameters
    ----------
    dvect1: data vector
        Data vector produced by :code:`2pt_pipeline`
    dvect2: data vector
        Data vector produced by :code:`WLpipe`
    nrow: int
        Number of rows for the pyplot.subplots
    ncol: int
        Number of columns for the pyplot.subplots
    nidx: list
        List of indexes for each subplot
    bins: list of tuples
        List of tuples (bin1,bin2) with the tomographic bin combinations
    ylabel: str
        The y-axis label
    legend: string
        Legend for subplots
    outname: str
        Output name for the plot
        
    Notes
    -----
    The output is a :code:`.png` file with the galaxy-galaxy comparison between
    :code:`2pt_pipeline` and :code:`WLpipe`.
    """
    plt.figure(1, figsize=(22,20))
    for i in range(len(bins)):
        angle1, value1, rbin = get_w_values(dvect1, bins[i])
        angle2, value2, _ = get_w_values(dvect2, bins[i])
        
        plt.subplot(nrow, ncol, nidx[i])     # two row, two columns, position 1
        if nidx[i]==6:
            plt.axis('off')
            plt.xlim(1e0, 450)
            plt.ylim(-0.5, 3.0)
            plt.text(2.5, 1.5, legend, fontsize=25)
            plt.text(2.5, 1.4, 'blue line=Y1 reference', color='C0', fontsize=25)
        else:
            plt.plot(angle1, value1, lw=2.5)    # y1
            plt.plot(angle2, value2, 'ko', mew=2, ms=12, mfc='None')
            plt.axhline(y=0, color='k', linestyle='--')
            plt.xlim(1e0, 450)
            plt.ylim(-0.5, 3.0)
            plt.xscale('log', nonposx='mask')
            plt.text(2.5, 2.5,'%d'%(rbin[0]), fontsize=25)
            plt.xlabel(r'$\theta$' + ' (arcmin)', fontsize=15)
            plt.ylabel(ylabel, fontsize=15) #(r'$\xi_{+}(\theta)$', fontsize=15)
            plt.tick_params(which='major', length=6, width=1, labelsize=15)
            plt.tick_params(which='minor', length=4, width=1)
    plt.savefig(outname, bbox_inches='tight', dpi=300)
    plt.close()


def plot_nz(nzs1, zs1, nzl1, zl1, nzs2, zs2, nzl2, zl2, nrow, ncol, outname):
    
    """
    Plot the comparison of n(z) of lens and sources between the
        :code:`2pt_pipeline` and :code:`WLpipe` on DES Y1 data.
        
    Parameters
    ----------
    nzs1: list of astropy.table
        The n(z) distribution of sources from :code:`2pt_pipeline`
    zs1: astropy.table
        Redshift of sources from :code:`2pt_pipeline`
    nzl1: list of astropy.table
        The n(z) distribution of lens from :code:`2pt_pipeline`
    zl1: astropy.table
        Redshift of lens from :code:`2pt_pipeline`

    nzs2: list of astropy.table
        The n(z) distribution of sources from :code:`WLpipe`
    zs2: astropy.table
        Redshift of sources from :code:`WLpipe`
    nzl2: list of astropy.table
        The n(z) distribution of lens from :code:`WLpipe`
    zl2: astropy.table
        Redshift of lens from :code:`WLpipe`
    nrow: int
        Number of rows for the pyplot.subplots
    ncol: int
        Number of columns for the pyplot.subplots
    outname: str
        Output name for the plot
    """
    plt.figure(1, figsize=(20,14))
    for i in range(nrow):
        plt.subplot(nrow, ncol, i+1)     # two row, two columns, position 1
        plt.xlabel('z', fontsize=20)
        plt.ylabel('n(z)', fontsize=20) #(r'$\xi_{+}(\theta)$', fontsize=15)
        plt.tick_params(which='both', length=6, width=1, labelsize=15)
        if i+1==1:
            plt.plot(zs2, nzs2[0], 'ko', mew=2, ms=5, mfc='None')
            plt.plot(zs2, nzs2[1], 'ko', mew=2, ms=5, mfc='None')
            plt.plot(zs2, nzs2[2], 'ko', mew=2, ms=5, mfc='None')
            plt.plot(zs2, nzs2[3], 'ko', mew=2, ms=5, mfc='None')
            plt.plot(zs1, nzs1[0], lw=2)    # y1
            plt.plot(zs1, nzs1[1], lw=2)
            plt.plot(zs1, nzs1[2], lw=2)
            plt.plot(zs1, nzs1[3], lw=2)
            plt.xlim(0, 2)
            plt.ylim(0, 5)
            plt.text(1.2, 4, 'points = WLpipe\nlines = Y1 reference', fontsize=15)
            plt.title('source redshifts', fontsize=15)

        elif i+1==2:
            plt.plot(zl2, nzl2[0], 'ko', mew=2, ms=5, mfc='None')
            plt.plot(zl2, nzl2[1], 'ko', mew=2, ms=5, mfc='None')
            plt.plot(zl2, nzl2[2], 'ko', mew=2, ms=5, mfc='None')
            plt.plot(zl2, nzl2[3], 'ko', mew=2, ms=5, mfc='None')
            plt.plot(zl2, nzl2[4], 'ko', mew=2, ms=5, mfc='None')
            plt.plot(zl1, nzl1[0], lw=2)    # y1
            plt.plot(zl1, nzl1[1], lw=2)
            plt.plot(zl1, nzl1[2], lw=2)
            plt.plot(zl1, nzl1[3], lw=2)
            plt.plot(zl1, nzl1[4], lw=2)
            plt.xlim(0, 1)
            plt.ylim(0, 10)
            plt.text(0.6, 8, 'points = WLpipe\nlines = Y1 reference', fontsize=15)
            plt.title('lens redshifts', fontsize=15)
    plt.savefig(outname, bbox_inches='tight', dpi=300)
    plt.close()

if __name__ == "__main__":

    # Path to the data vectors
    #wlpipe_outpath='/data/des61.b/data/mariaeli/demo_scott/output/'
    #desy1_outpath = '/data/des61.b/data/mwang/data/DES/Y1/data_vectors/'
    wlpipe_outpath='/Users/maria/current-work/wlpipe/'
    desy1_outpath = '/Users/maria/current-work/wlpipe/'

    # Plotting xip, bin1==bin2
    y1_xip_dvect, xip_bin1, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 2)
    wl_xip_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 2)
    xi_nrow = len(xip_bin1)
    xi_nidx = [1,2,3,4,5,6,7,9,10,11,13]
    xi_binning = [(0,3), (1,3), (2,3), (3,3), (0,2), (1,2), (2,2), (0,1), (1,1), (0,0), (0,0)]
    
    plot_xi(y1_xip_dvect, wl_xip_dvect, xip_bin1, xi_nrow, xi_nidx, xi_binning, r'$\xi_{+}(\theta)$',
            r'shear-shear $\xi_{+}$'+'\n\n'+'points=WLpipe'+'\n', 'desy1_xip.png')

    # Plotting xim, bin1==bin2
    y1_xim_dvect, xim_bin1, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 3)
    wl_xim_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 3)
    
    plot_xi(y1_xim_dvect, wl_xim_dvect, xim_bin1, xi_nrow, xi_nidx, xi_binning, r'$\xi_{-}(\theta)$',
            r'shear-shear $\xi_{-}$'+'\n\n'+'points=WLpipe'+'\n', 'desy1_xim.png')

    # Plotting gammat, bin1!=bin2
    y1_gammat_dvect, gammat_bin1, gammat_bin2 = get_fits(desy1_outpath, '2pt_NG_1101.fits', 4)
    wl_gammat_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 4)
    gammat_ncol = len(gammat_bin1)
    gammat_nrow = len(gammat_bin2)
    gammat_nidx = range(1,21)
    gammat_binning = []
    for i in range(1,5):
        for j in range(1,6):
            gammat_binning.append((j,i))

    plot_gammat(y1_gammat_dvect, wl_gammat_dvect, gammat_nrow, gammat_ncol, gammat_nidx, gammat_binning,
                r'$\theta \, \gamma_{t}\,$'+r'($10^{-2}$ arcmin)', 'desy1_gammat.png')

    # Plotting w, bin1==bin2
    y1_w_dvect, w_bin1, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 5)
    wl_w_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 5)
    w_nrow = 2
    w_ncol = 3
    w_binning = [(i,i) for i in range(1,6)] + [(5,5)]
    w_nidx = range(1,7)

    plot_w(y1_w_dvect, wl_w_dvect, w_nrow, w_ncol, w_nidx, w_binning, r'$\theta \,w \,$'+r'(arcmin)',
           r'galaxy-galaxy $w(\theta)$'+'\n\n'+'points=WLpipe'+'\n', 'desy1_w.png')

    # Getting nz_source, bin1==bin2
    y1_nzsource_dvect, nzsource_bin1, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 6)
    wl_nzsource_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 6)

    # Getting nz_lens, bin1==bin2
    y1_nzlens_dvect, nzlens_bin1, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 7)
    wl_nzlens_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 7)

    # Plotting n(z)
    y1_zs1, y1_zs2, y1_zs3, y1_bs1, y1_bs2, y1_bs3, y1_bs4, y1_bs5 = get_nz_values(y1_nzsource_dvect)
    y1_zl1, y1_zl2, y1_zl3, y1_bl1, y1_bl2, y1_bl3, y1_bl4, y1_bl5 = get_nz_values(y1_nzlens_dvect)

    wl_zs1, wl_zs2, wl_zs3, wl_bs1, wl_bs2, wl_bs3, wl_bs4, wl_bs5 = get_nz_values(wl_nzsource_dvect)
    wl_zl1, wl_zl2, wl_zl3, wl_bl1, wl_bl2, wl_bl3, wl_bl4, wl_bl5 = get_nz_values(wl_nzlens_dvect)

    plot_nz([y1_bs1, y1_bs2, y1_bs3, y1_bs4], y1_zs2, [y1_bl1, y1_bl2, y1_bl3, y1_bl4, y1_bl5], y1_zl2,
            [wl_bs1, wl_bs2, wl_bs3, wl_bs4], wl_zs2, [wl_bl1, wl_bl2, wl_bl3, wl_bl4, wl_bl5], wl_zl2,
            2, 1, 'desy1_dndz.png')



