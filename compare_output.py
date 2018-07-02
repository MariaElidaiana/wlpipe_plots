#!/usr/env python
"""
.. module:: compare_output
:synopsis: Script to check the data vector results from WLpipe and compare with DES Y1.
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

def dvector_comparison(dvectx1, dvectx2, dvecty1, dvecty2, dec1, dec2, errmens1, errmens2, dvname):
    """
    Reads the data vector as astropy.table. Setup the precision
    and error mensagens for the :code:`np.testing.assert_array_almost_equal`
    function. Test if the two data-vectors have the same x,y-axis values.
        
    Parameters
    ----------
    dvectx1: astropy.table
        x-axis of the data vector 1
    dvectx2: astropy.table
        x-axis of the data vector 2
    dvecty1: astropy.table
        y-axis of the data vector 1
    dvecty2: astropy.table
        y-axis of the data vector 2
    dec1: int
        Desired precision for x-axis of the data vectors
    dec2: int
        Desired precision for y-axis of the data vectors
    errmens1: str
        Error message if x-axis data vectors comparison fail
    errmens2: str
        Error message if y-axis data vectors comparison fail
    dvname: str
        Name of the data vector quantity that is being compared
        
    Returns
    -------
        
    Notes
    -----
    Change the decimal precision to check when the data vectors values
    comparison fail.
    """

    try:
        np.testing.assert_array_almost_equal(np.array(dvectx1), np.array(dvectx2), decimal=dec1, err_msg=errmens1)
        t1=True
    except AssertionError as e:
        print e
        print '-'*80

    try:
        np.testing.assert_array_almost_equal(np.array(dvecty1), np.array(dvecty2), decimal=dec2,  err_msg=errmens2)
        t2=True
    except AssertionError as e:
        print e
        print '='*80

    try:
        if t1&t2:    print '- Checking: ' +dvname+' data vectors are the same! (x,y)-axis precision: ('+str(dec1)+','+str(dec2)+') decimals.'
        print '='*80
    except UnboundLocalError:
        pass


if __name__ == "__main__":

    # Path to the data vectors
    #wlpipe_outpath='/data/des61.b/data/mariaeli/demo_scott/output/'
    #desy1_outpath = '/data/des61.b/data/mwang/data/DES/Y1/data_vectors/'
    wlpipe_outpath='/Users/maria/current-work/wlpipe/'
    desy1_outpath = '/Users/maria/current-work/wlpipe/'

    # xi_plus checking
    y1_xip_dvect, _, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 2)
    wl_xip_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 2)

    dvector_comparison(y1_xip_dvect['ANG'], wl_xip_dvect['ANG'],
                       y1_xip_dvect['VALUE'], wl_xip_dvect['VALUE'],
                       2, 6,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): xi_plus angles are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): xi_plus are not the same:\n',
                       'xi_plus')

    # xi_minus checking
    y1_xim_dvect, _, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 3)
    wl_xim_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 3)

    dvector_comparison(y1_xim_dvect['ANG'], wl_xim_dvect['ANG'],
                       y1_xim_dvect['VALUE'], wl_xim_dvect['VALUE'],
                       2, 6,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): xi_minus angles are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): xi_minus are not the same:\n',
                       'xi_minus')

    # gamma_t checking
    y1_gammat_dvect, _, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 4)
    wl_gammat_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 4)
    dvector_comparison(y1_gammat_dvect['ANG'], wl_gammat_dvect['ANG'],
                       y1_gammat_dvect['VALUE'], wl_gammat_dvect['VALUE'],
                       2, 5,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): gamma_t angles are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): gamma_t are not the same:\n',
                       'gamma_t')

    # w checking
    y1_w_dvect, _, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 5)
    wl_w_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 5)
    dvector_comparison(y1_w_dvect['ANG'], wl_w_dvect['ANG'],
                       y1_w_dvect['VALUE'], wl_w_dvect['VALUE'],
                       2, 4,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): w angles are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): w are not the same:\n',
                       'w')

    # nzsource checking
    y1_nzsource_dvect, _, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 6)
    wl_nzsource_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 6)

    dvector_comparison(y1_nzsource_dvect['Z_MID'], wl_nzsource_dvect['Z_MID'],
                       y1_nzsource_dvect['BIN1'], wl_nzsource_dvect['BIN1'],
                       4, 2,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_source_bin1 z_mid are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_source_bin1 are not the same:\n',
                       'nz_source_bin1')

    dvector_comparison(y1_nzsource_dvect['Z_MID'], wl_nzsource_dvect['Z_MID'],
                       y1_nzsource_dvect['BIN2'], wl_nzsource_dvect['BIN2'],
                       4, 2,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_source_bin2 z_mid are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_source_bin2 are not the same:\n',
                       'nz_source_bin2')

    dvector_comparison(y1_nzsource_dvect['Z_MID'], wl_nzsource_dvect['Z_MID'],
                       y1_nzsource_dvect['BIN3'], wl_nzsource_dvect['BIN3'],
                       4, 2,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_source_bin3 z_mid are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_source_bin3 are not the same:\n',
                       'nz_source_bin3')

    dvector_comparison(y1_nzsource_dvect['Z_MID'], wl_nzsource_dvect['Z_MID'],
                       y1_nzsource_dvect['BIN4'], wl_nzsource_dvect['BIN4'],
                       4, 2,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_source_bin4 z_mid are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_source_bin4 are not the same:\n',
                       'nz_source_bin4')

    # nzlens checking
    y1_nzlens_dvect, _, _ = get_fits(desy1_outpath, '2pt_NG_1101.fits', 7)
    wl_nzlens_dvect, _, _ = get_fits(wlpipe_outpath, 'desy1_dvect.fits', 7)

    dvector_comparison(y1_nzlens_dvect['Z_MID'], wl_nzlens_dvect['Z_MID'],
                       y1_nzlens_dvect['BIN1'], wl_nzlens_dvect['BIN1'],
                       4, 1,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin1 z_mid are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin1 are not the same:\n',
                       'nz_lens_bin1')

    dvector_comparison(y1_nzlens_dvect['Z_MID'], wl_nzlens_dvect['Z_MID'],
                       y1_nzlens_dvect['BIN2'], wl_nzlens_dvect['BIN2'],
                       4, 1,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin2 z_mid are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin2 are not the same:\n',
                       'nz_lens_bin2')

    dvector_comparison(y1_nzlens_dvect['Z_MID'], wl_nzlens_dvect['Z_MID'],
                       y1_nzlens_dvect['BIN3'], wl_nzlens_dvect['BIN3'],
                       4, 1,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin3 z_mid are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin3 are not the same:\n',
                       'nz_lens_bin3')

    dvector_comparison(y1_nzlens_dvect['Z_MID'], wl_nzlens_dvect['Z_MID'],
                       y1_nzlens_dvect['BIN4'], wl_nzlens_dvect['BIN4'],
                       4, 1,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin4 z_mid are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin4 are not the same:\n',
                       'nz_lens_bin4')

    dvector_comparison(y1_nzlens_dvect['Z_MID'], wl_nzlens_dvect['Z_MID'],
                       y1_nzlens_dvect['BIN5'], wl_nzlens_dvect['BIN5'],
                       4, 1,
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin5 z_mid are not the same:\n',
                       '\nY1_2pt_pipeline (x) vs Y1_WLpipe (y): nz_lens_bin5 are not the same:\n',
                       'nz_lens_bin5')




