import numpy as np
import pint
units = pint.UnitRegistry()
import astropy.io.fits as fits
import basis
import matplotlib.pyplot as plt
from datetime import datetime

# Script that reads in a set of CodeV Pupil map OPD files across fields and instruments and fits Zernikes to the OPD
# distribution at each field point and then Legendres to the variation of each Zernike coefficient across field.  The
# Resulting table of Legendres coefficients for each Zernike term is written to a .fits table which is then used in
# WebbPSF to model field dependence.

def read_codeV_data_file(filename):
    headerlines = 16
    with open(filename) as fh:
        for i, line in enumerate(fh):
            if i is headerlines: break
            if i is 6:
                x_field_deg = float(str.split(line)[3])
            if i is 7:
                y_field_deg = float(str.split(line)[3])
    data = np.loadtxt(filename, skiprows=headerlines)
    data_defined = data != -99999

    return data, data_defined, (x_field_deg, y_field_deg)


def main():
    num_fields_x = 21
    num_fields_y = 21

    num_wf_pts_x = 256
    num_wf_pts_y = 256

    instruments = [
        'fgs',
        'nircam',
        'miri',
        'nirspec',
        'niriss'
    ]

    # Zernike fitting order for each instrument
    order= {'fgs': 15,
            'nircam': 15,
            'miri': 15,
            'nirspec': 15,
            'niriss': 15
            }

    # Define the origin of the CodeV coordinate system in V2/V3 space
    v2_origin_degrees = 0
    v3_origin_degrees = -468/3600
    # How are the signs of the coordinates related between V2/V3 and CodeV
    v2_sign = 1
    v3_sign = -1

    wavelength_nm = 2120

    # Fitting order of Legendres for fitting field variation
    legendre_order = 15
    legendres = basis.Legendre2DBasis(legendre_order, num_fields_y, num_fields_x)

    # Read in the wavefront at at the on-axis field point
    oa_data, oa_data_defined, oa_field = read_codeV_data_file('wavefront_on_axis.dat')
    oa_data = oa_data * wavelength_nm

    # Loops over the various instruments
    for instrument in instruments:
        try:
            if zernikes.order != order[instrument]:
                zernikes = basis.NollZernikeBasis(order[instrument], num_wf_pts_x, num_wf_pts_y)
                num_coeffs = zernikes.numpolys
        except NameError:
            zernikes = basis.NollZernikeBasis(order[instrument], num_wf_pts_x, num_wf_pts_y)
            num_coeffs = zernikes.numpolys

        # Read in wavefront at reference point
        filename = f'{instrument}/wavefront_ref.dat'
        ref_data, ref_data_defined, ref_field = read_codeV_data_file(filename)
        ref_data = ref_data * wavelength_nm  # CodeV gives wavefront in waves

        # Initialize some arrays
        fields = np.zeros((num_fields_x, num_fields_y, 2))
        data_defined_arr =  np.zeros((num_fields_x, num_fields_y, num_wf_pts_x, num_wf_pts_y)).astype(bool)

        coeffs_local = np.zeros((num_fields_y, num_fields_x, num_coeffs))
        coeffs_global = np.zeros((num_fields_y, num_fields_x, num_coeffs))

        print(instrument)

        # Loop over the files for all fields doing Zernike fits and storing the Zernike data
        for fx_index in range(0, num_fields_x):

            for fy_index in range(0, num_fields_y):
                filename = f'{instrument}/wavefront_fx{fx_index + 1}_fy{fy_index + 1}.dat'
                print(filename)
                # Read in the CodeV wavefront data, OPD in waves in CodeV
                cur_data, cur_data_defined, cur_field_pt = read_codeV_data_file(filename)
                cur_data = cur_data * wavelength_nm

                # Keep track of where data is good for each wavefront.  And together the current field, reference field
                # and on-axis field so we only consider data that is valid in all of them.
                data_defined_arr[fx_index, fy_index, :, :] = np.logical_and(cur_data_defined, ref_data_defined)
                data_defined_arr[fx_index, fy_index, :, :] = np.logical_and(data_defined_arr[fx_index, fy_index, :, :], oa_data_defined)

                # Store the location of the current field point
                fields[fx_index, fy_index, 0] = cur_field_pt[0]
                fields[fx_index, fy_index, 1] = cur_field_pt[1]

                # Subtract off local reference (~center of current instrument) and global (0,0 CodeV field point)
                # reference wavefronts from current wavefront
                cur_data_local = cur_data - ref_data
                cur_data_global = cur_data - oa_data

                # Define wavefront objects using the CodeV wavefront data
                cv_wavefront_global = basis.PointByPointWavefront(cur_data_global,
                                                           data_defined_arr[fx_index, fy_index, :, :], 2, 2120,
                                                           basis.units.nm, wf_basis=zernikes)
                cv_wavefront_local = basis.PointByPointWavefront(cur_data_local,
                                                           data_defined_arr[fx_index, fy_index, :, :], 2, 2120,
                                                           basis.units.nm, wf_basis=zernikes)

                # Get the coefficients for the wavefront from their wavefront objects and set piston, tip and tilt
                # contributions to zero.
                coeffs_global[fy_index, fx_index, :] = cv_wavefront_global.coeffs

                # ttp_global = np.zeros(zernikes.numpolys)
                # ttp_global[0:3] = coeffs_global[fy_index, fx_index, 0:3]
                coeffs_global[fy_index, fx_index, 0:3] = 0
                coeffs_local[fy_index, fx_index, :] = cv_wavefront_local.coeffs
                coeffs_local[fy_index, fx_index, 0:3] = 0


        # Fit the field dependent coefficients to Legendres.  Each Zernike varying across field gets its own Legendre fit

        #Initialize column lists used to store data in fits file
        fits_column_list_local = []
        fits_column_list_global = []
        #Initialize arrays for the current set of coefficients
        legendre_coeffs_local = np.zeros((legendres.numpolys, num_coeffs))
        legendre_coeffs_global = np.zeros((legendres.numpolys, num_coeffs))

        #Loop over each Zernike, doing the fits, etc
        for index in range(0, num_coeffs):
            # Get array of current Zernike coefficient across field
            cur_coeffs_local = coeffs_local[:,:, index]
            cur_coeffs_global = coeffs_global[:, :, index]
            # Set the Zernike coefficient variation up as if it was a wavefront.  Of course the array here is varying
            # across the field plane and not the pupil plane, so we're using this in a non-standard way.  Also,
            # using the legendre polynomials at the basis set.
            legendre_var_local = basis.PointByPointWavefront(cur_coeffs_local, np.ones_like(cur_coeffs_local).astype(bool), 4,
                                                             wavelength_nm, basis.units.nm, wf_basis=legendres)
            legendre_var_global = basis.PointByPointWavefront(cur_coeffs_global, np.ones_like(cur_coeffs_global).astype(bool), 4,
                                                             wavelength_nm, basis.units.nm, wf_basis=legendres)
            # Pull out the Legendre coefficients from the wavefront object
            legendre_coeffs_local[:, index] = legendre_var_local.coeffs
            legendre_coeffs_global[:, index] = legendre_var_global.coeffs

            # Store the Legendre coefficients for the current Zernike term as a fits column and store in column list.
            cur_fits_col_local = fits.Column(name=f'Zernike {index} Legendres', format='D', array=legendre_coeffs_local[:,index])
            fits_column_list_local.append(cur_fits_col_local)
            cur_fits_col_global = fits.Column(name=f'Zernike {index} Legendres', format='D', array=legendre_coeffs_global[:,index])
            fits_column_list_global.append(cur_fits_col_global)

        # Set things up to write coefficient data to .fits file

        # Determine the limits of the field region over which we the current model is valid
        min_xfield = np.min(fields[:, :, 0])
        max_xfield = np.max(fields[:, :, 0])
        min_yfield = np.min(fields[:, :, 1])
        max_yfield = np.max(fields[:, :, 1])

        # Instantiate a fits object and populate header fields
        hdu_primary = fits.PrimaryHDU([])
        header = hdu_primary.header
        header['instr'] = instrument
        header.comments['instr'] = 'Name of the instrument field described'
        header['minxfie'] = min_xfield * 60
        header.comments['minxfie'] = 'Min valid X field coordinate'
        header['maxxfie'] = max_xfield * 60
        header.comments['maxxfie'] = 'Max valid X field coordinate'
        header['minyfie'] = min_yfield * 60
        header.comments['minyfie'] = 'Min valid Y field coordinate'
        header['maxyfie'] = max_yfield * 60
        header.comments['maxyfie'] = 'Max valid Y field coordinate'
        header['wfbasis'] = 'Noll Zernikes'
        header.comments['wfbasis'] = 'Basis polynomial set used for OPD in the pupil'
        header['fiebasis'] = 'Legendre Polynomials'
        header.comments['fiebasis'] = 'Basis set for variation of pupil coeffs with field'
        header['ncoefwf'] = zernikes.numpolys
        header.comments['ncoefwf'] = 'No. of polys for the pupil OPD distribution'
        header['ncoeffie'] = legendres.numpolys
        header.comments['ncoeffie'] = 'No. of polys for variation of pupil coeffs with field'
        header['fieorder'] = legendres.order
        header.comments['fieorder'] = 'Max. order of polynomials for field variation'
        header['v2origin'] = v2_origin_degrees * 60
        header.comments['v2origin'] = 'Location of origin of V2 coordinate system'
        header['v3origin'] = v3_origin_degrees * 60
        header.comments['v3origin'] = 'Location of origin of V3 coordinate system'
        header['v2sign'] = v2_sign
        header.comments['v2sign'] = 'Sign of the coord. sys. used in this file wrt to V2'
        header['v3sign'] = v3_sign
        header.comments['v3sign'] = 'Sign of the coord. sys. used in this file wrt to V3'
        header['refptx'] = ref_field[0] * 60
        header.comments['refptx'] = 'X coordinate of reference point (center) of ins. field'
        header['refpty'] = ref_field[1] * 60
        header.comments['refpty'] = 'Y coordinate of reference point (center) of ins. field'
        header['fangunit'] = 'arcmin'
        header.comments['fangunit'] = '(Angular) units of field coordinates'
        header['opdunit'] = 'nm'
        header.comments['opdunit'] = 'Units used for the OPDs'
        now = datetime.now()
        header['created'] = now.strftime('%m/%d/%Y, %H:%M:%S')
        header.comments['created'] = 'File creation date and time'

        # Put local and global fitting data in to fits hdu as tables
        col_defs_local = fits.ColDefs(fits_column_list_local)
        col_defs_global = fits.ColDefs(fits_column_list_global)
        hdu_table_local = fits.BinTableHDU.from_columns(col_defs_local)
        hdu_table_global = fits.BinTableHDU.from_columns(col_defs_global)

        # Populate the hdu list of our fits object and write the data out.
        hdu_list = fits.HDUList([hdu_primary, hdu_table_local, hdu_table_global])
        hdu_list.writeto(f'{instrument}/field_dep_table_{instrument}.fits', overwrite='True')

if __name__ == '__main__':
    main()