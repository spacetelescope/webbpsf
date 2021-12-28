import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
import math

# This is the start of an experimental prototype object oriented code for representing wavefronts.  The intent is to
# represent a wavefront in such a way that calling code doesn't need to care about the representation and can just ask
# for the data as an array of OPD values on a regular grid, polynomial coefficients, or other without caring about
# the internal representation or worrying about fitting or expanding polynomials out, etc.
# This code mostly works as is, but there are some features not fully implemented or implemented stupidly.  No effort
# has been made to optimize it, so it is crazy slow.  It also tries to implement unit tracking and conversions for OPD
# and phase values, but this doesn't work quite right either.  However, the polynomial representations and fitting
# parts seem to work well since they are ported from older code.
# Use this at your own risk (or better, don't.)  When I get time perhaps it will be get better, but since there are
# other options most people should probably use those.
# This was mostly written as an exercise and for my own use to explore some ideas I had for representing things
# cleanly.  Greg Brady, January 2021.

# Disable Pint's old fallback behavior (must come before importing Pint)
import os

os.environ['PINT_ARRAY_PROTOCOL_FALLBACK'] = "0"

import pint

units = pint.UnitRegistry()
Q_ = units.Quantity

# Silence NEP 18 warning
import warnings

with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    Q_([])

import copy


def embed(n, m):
    # Return vector of indices needed to pull an array of length m from the center of an array of (larger) length n, or
    # conversely insert an array of length m into an array of length n.
    # Idea:  there should be a way to make a class extending standard array class to handle this automatically.
    if n < m:
        raise ValueError("Argument n is less than argument m")
    emvec = np.fix(n / 2) + range(-np.fix(m / 2).astype('int64'), np.fix((m - 1) / 2).astype('int64') + 1)
    return emvec.astype('int')


class BaseBasis:
    """Base class for basis functions that are used to represent phase distributions of wavefront polynomials"""

    # def __init__(self, pts_x, pts_y):
    def __init__(self, *args, **kwargs):

        """Base class constructor that sets attributes that every basis set will need and creates
       a 2D grid of points on which the basis set will be represented.  It is assumes that the
       basis set is scaled so that the array will generally run from -1 to 1 in x and y.  In
       cases where the number of points is odd, the arrays will be centered around 0 and run from
       -1 to 1 exactly.  In cases where the number of points is even, the arrays will be centered
       around (pts_x/2, pts_y/2) and the first array indices (i.e. [0]) will have values one step
       less than -1.
       Arguments:
           pts_x:  number of points representing basis in x dimension
           pts_y:  number of points representing basis in y dimension"""
        self._pts_x = args[0]
        self._pts_y = args[1]

        # Calculate positive and negative maximum number of points to have
        min_x = -np.fix(self._pts_x / 2)
        max_x = np.fix((self._pts_x + 1) / 2)
        min_y = -np.fix(self._pts_y / 2)
        max_y = np.fix((self._pts_y + 1) / 2)

        # Step sizes to get -1 to 1 range
        self._step_x = 1 / np.fix((self._pts_x - 1) / 2)
        self._step_y = 1 / np.fix((self._pts_y - 1) / 2)

        # Create open arrays for coordinates and scales
        self._ynorm, self._xnorm = np.ogrid[min_y:max_y, min_x:max_x]
        self._xnorm = self._xnorm * self._step_x
        self._ynorm = self._ynorm * self._step_y

        self._extent = np.array(
            [min_x * self._step_x, (max_x - 1) * self._step_x, min_y * self._step_y, (max_y - 1) * self._step_y])

        # Create placeholder array that represents where data is defined (i.e. what will generally
        # be used to represent an aperture or the like.
        self._data_defined = np.ones((self._pts_y, self._pts_x))
        # Typically the first row and column won't have data defined if there are an even number of points in the array
        if (self._pts_x % 2) == 0:
            self._data_defined[:, 0] = 0
        if (self._pts_y % 2) == 0:
            self._data_defined[0, :] = 0

    @property
    def step_x(self):
        return self._step_x

    @property
    def step_y(self):
        return self._step_y

    @property
    def xnorm(self):
        return self._xnorm

    @property
    def ynorm(self):
        return self._ynorm

    @property
    def data_defined(self):
        return self._data_defined

    @property
    def extent(self):
        return self._extent

    @property
    def pts_x(self):
        return self._pts_x

    @pts_x.setter
    def pts_x(self, new_pts_x):
        self.__init__(new_pts_x, self._pts_y)

    @property
    def pts_y(self):
        return self._pts_y

    @pts_y.setter
    def pts_y(self, new_pts_y):
        self.__init__(self._pts_x, new_pts_y)


class PolynomialBasis(BaseBasis):
    """Class that implements polynomial basis sets"""

    def __init__(self, *args, **kwargs):

        self._order = args[0]
        BaseBasis.__init__(self, args[1], args[2])
        # Set up dummy polynomial array and number of polynomials
        self._polynomials = np.ones((self._order, args[1], args[2]))
        self._numpolys = args[0]
        self._rms_factors = np.ones(self._numpolys)
        self._term_order = None

        self._normalization = kwargs.pop('normalization', 'natural')

    # @property
    # def pts_x(self):
    #     return self._pts_x

    @property
    def order(self):
        return self._order

    @order.setter
    def order(self, new_order):
        self.__init__(new_order, self.pts_x, self.pts_y)

    @property
    def pts_y(self):
        return super().pts_y

    @pts_y.setter
    def pts_y(self, new_pts_y):
        self.__init__(self._order, self._pts_x, new_pts_y)

    @property
    def pts_x(self):
        return super().pts_x

    @pts_x.setter
    def pts_x(self, new_pts_x):
        self.__init__(self._order, new_pts_x, self._pts_y)

    @property
    def rms_factors(self):
        return self._rms_factors

    @property
    def polynomials(self):
        if self._normalization == 'natural':
            return self._polynomials
        elif self._normalization == 'unit_pv':
            raise NotImplementedError(
                'Unit PV normalization is not implemented.  Previous implementation broke orthogonality')
            # TODO:  FIX THIS!!!  An offset like this to force things to run from -1 to 1 will break orthogonality, and
            #  must not be done.  Think about what to do instead.  Check Noll paper.  Possible just make PV truly unity
            #  and don't care about what the absolute max and min values are.  Also implement property decorator and
            #  change decorator for changing the output normalization.
            # scale = [1]
            # offset = [0]
            # for polynomial in self._polynomials[1:, :, :]:
            #     max_val = np.max(polynomial[self.data_defined])
            #     min_val = np.min(polynomial[self.data_defined])
            #     cur_scale = 2 / (max_val - min_val)
            #     scale.append(cur_scale)
            #     offset.append(-(min_val * cur_scale + 1))
            # offset = np.array(offset)
            #

            # return (np.einsum('ijk, i->ijk', self._polynomials, scale)) + offset.reshape(self._numpolys, 1, 1)
        elif self._normalization == 'unit_rms':
            return np.einsum('ijk, i->ijk', self._polynomials, 1 / self._rms_factors)
        else:
            raise ValueError('Invalid normalization selected')

    @property
    def numpolys(self):
        return self._numpolys

    @numpolys.setter
    def numpolys(self, new_numpolys):
        polynomials_copy = self._polynomials
        self._polynomials = np.zeros((new_numpolys, self._pts_x, self._pts_y))
        self._polynomials[0:self._numpolys, :, :] = polynomials_copy
        self._numpolys = new_numpolys

    def project(self, in_data):
        coeffs = np.zeros(self.numpolys)
        for polynomial, index in zip(self._polynomials, range(0, self.numpolys)):
            coeff = np.einsum('i,i', polynomial[self.data_defined], in_data[self.data_defined]) / \
                    np.einsum('i,i->', polynomial[self.data_defined], polynomial[self.data_defined])
            coeffs[index] = coeff

        return coeffs

    def polyfit(self, in_data, in_data_defined):
        gs_poly = GramSchmidtBasis(in_data_defined, self)
        # for poly in gs_poly.polynomials:
        #     fig, ax = plt.subplots()
        #     im = plt.imshow(poly * gs_poly.data_defined)
        #     fig.colorbar(im, ax=ax)
        #     plt.show()
        gs_coeffs = gs_poly.project(in_data)

        # Calculate terms to go between original basis and new orthogonal basis
        c = np.zeros([self.numpolys, self.numpolys])
        for r in range(0, self.numpolys):
            c[r, r] = 1
            for index in range(0, r):
                c[r, index] = 0
                for s in range(0, r - index):
                    c[r, index] = c[r, index] + gs_poly._d[r, r - s - 1] * c[r - s - 1, index]

        # Calculate the coefficients for the original basis from the coefficients for the orthogonal basis
        fit_coeffs = np.einsum('ij,i->j', c, gs_coeffs)  # Multiply coefficients by conversion matrix c

        # Loop based implementation, left here so I can understand what's going on above better.
        # fit_coeffs = np.zeros(self.numpolys)
        # for r in range(0, self.numpolys):
        #     fit_coeffs[r] = gs_coeffs[r]
        #     for index in range(r + 1, self.numpolys):
        #         fit_coeffs[r] = fit_coeffs[r] + gs_coeffs[index] * c[index, r]
        # fit_coeffs[self.numpolys - 1] = gs_coeffs[self.numpolys - 1]
        return fit_coeffs

    def show_polynomials(self, plot_terms=None, individual=False):
        if not plot_terms:
            plot_terms = range(0, self._numpolys)

        if individual:
            for polynomial in self._polynomials[plot_terms, :, :]:
                fig, ax = plt.subplots()
                im = ax.imshow(polynomial * self.data_defined, extent=self._extent)
                fig.colorbar(im, ax=ax)
        else:
            min_val = np.min(self._polynomials[:, np.transpose(self._data_defined)])
            max_val = np.max(self._polynomials[:, np.transpose(self._data_defined)])

            max_order = int(np.max(self._term_order[plot_terms]))
            plotted = np.zeros((max_order + 1, max_order + 1), dtype=bool)
            fig, ax = plt.subplots(max_order + 1, max_order + 1)
            for polynomial, index in zip(self.polynomials[plot_terms, :, :], plot_terms):
                col = index - int(self._term_order[index] * (self._term_order[index] + 1) / 2)
                im = ax[self._term_order[index], col].imshow(polynomial * np.transpose(self.data_defined), extent=self._extent,
                                                             vmin=min_val, vmax=max_val)
                ax[self._term_order[index], col].axis('off')
                plotted[self._term_order[index], col] = True

            for cur_ax in ax[np.nonzero(np.logical_not(plotted))]:
                cur_ax.axis('off')
            fig.subplots_adjust(right=0.825)
            cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
            fig.colorbar(im, cax=cbar_ax)


class ZernikeBasis(PolynomialBasis):
    """Class that implements Zernike polynomial Basis"""

    def __urect__(self):
        """Calculates Zernike polynomial to an arbitrary order on a rectangular grid"""

        # Calculate zernike polynomials up to specified order in rectangular coordinates.  Follows Malacara 2nd Ed.,
        # pp461 - 472

        x_powers = np.zeros((self._order + 1, self._ynorm.shape[0], self._xnorm.shape[1]))
        y_powers = np.zeros((self._order + 1, self._ynorm.shape[0], self._xnorm.shape[1]))

        # Calculate powers of x up to the order of the polynomial required
        x_powers[0, :, :] = 1
        y_powers[0, :, :] = 1
        x_powers[1, :, :] = self._xnorm
        y_powers[1, :, :] = self._ynorm

        for index in range(2, self.order + 1):
            x_powers[index, :, :] = np.einsum('..., ...', x_powers[index - 1, :, :], x_powers[1, :, :])
            y_powers[index, :, :] = np.einsum('..., ...', y_powers[index - 1, :, :], y_powers[1, :, :])

        # Calculate number of polynomials up to the given order(triangular number)
        self._numpolys = int((self._order + 1) * (self._order + 2) / 2)

        self._polynomials = np.zeros((self._numpolys, self._ynorm.shape[0], self._xnorm.shape[1]))

        # array to keep track of which terms are rotationally invariant
        self._radial_terms = np.zeros(self.numpolys, dtype=bool)

        # array to keep track of the order of the polynomial
        self._term_order = np.zeros(self._numpolys).astype(int)

        # array to keep track of factors to renormalize terms for unit RMS WF error
        self._rms_factors = np.zeros(self._numpolys)

        # loop over the polynomials
        for r_index in range(0, self._numpolys):
            # Calculate the order of the current polynomial
            n = int(np.ceil((-3 + (9 + 8 * r_index) ** 0.5) / 2))

            self._term_order[r_index] = n

            m = r_index - (n + 1) * n // 2

            if (n - 2 * m) > 0:  # Sine function
                p = 1
                j_sumlim = m
                if (n % 2) == 0:
                    q = ((n - 2 * m) // 2 - 1)  # n even
                else:
                    q = (n - 2 * m - 1) // 2  # n odd
            else:  # Cos function
                p = 0
                j_sumlim = n - m
                if (n % 2) == 0:
                    q = (2 * m - n) // 2  # n even
                else:
                    q = (2 * m - n - 1) // 2  # n odd

            # Initialize the current Zernike
            self._polynomials[r_index, :, :] = np.zeros((self._ynorm.shape[0], self._xnorm.shape[1]))

            # Loop over sum indices
            for i_index in range(0, q + 1):
                for j_index in range(0, j_sumlim + 1):
                    for k_index in range(0, j_sumlim - j_index + 1):
                        self._polynomials[r_index, :, :] = self._polynomials[r_index, :, :] \
                                                           + (-1) ** (i_index + j_index) \
                                                           * sp.comb(abs(2 * m - n), 2 * i_index + p) \
                                                           * sp.comb(j_sumlim - j_index, k_index) \
                                                           * (sp.gamma(n - j_index + 1)
                                                              / (sp.gamma(j_index + 1) * sp.gamma(m - j_index + 1)
                                                                 * sp.gamma(n - m - j_index + 1))) \
                                                           * x_powers[n - 2 * (i_index + j_index + k_index) - p, :, :] \
                                                           * y_powers[2 * (i_index + k_index) + p, :, :]

                        # Calculate the multiplicative factor to normalize the Zernike for unit RMS error
                        # useful for calculating RMS wavefront error.
                        if n - 2 * m == 0:
                            # self._rms_factors[r_index] = 1 / np.sqrt((n + 1))
                            self._rms_factors[r_index] = 1 / np.sqrt((n + 1))

                        else:
                            # self._rms_factors[r_index] = 1 / np.sqrt(2 * (n + 1))
                            self._rms_factors[r_index] = 1 / np.sqrt(2 * (n + 1))

                        # Flag radial terms
                        if n - 2 * m == 0:
                            self._radial_terms[r_index] = True

    def __init__(self, *args, **kwargs):

        input_norm = kwargs.pop('normalization', 'natural')
        if input_norm not in {'natural', 'unit_pv', 'unit_rms'}:
            raise ValueError('Selected normalization not supported')
        PolynomialBasis.__init__(self, *args, normalization=input_norm, **kwargs)

        # Set the data defined mask as a unit circle
        self._r = np.sqrt(np.power(self._xnorm, 2) + np.power(self._ynorm, 2))
        self._data_defined = self._r <= 1

        # Calculate the Zernikes
        self.__urect__()

    @property
    def r(self):
        return self._r

    # @property
    # def numpolys(self):
    #     return self._numpolys

    @property
    def radial_terms(self):
        return self._radial_terms

    @property
    def term_order(self):
        return self._term_order

class NollZernikeBasis(ZernikeBasis):
    """Class that calculates the Noll ordering/normalization of the Zernikes"""

    def __calc_noll_zerns__(self):

        j_index = 0

        r = np.sqrt(self._xnorm ** 2 + self._ynorm ** 2)
        theta = np.arctan2(self._ynorm, self._xnorm)

        cos_sin = False

        for n in range(0, self.order + 1):
            for m in range((n % 2), (n + 1), 2):
                r_n_m = np.zeros_like(r)
                for s in range(0, int((n - m) / 2) + 1):
                    r_n_m = r_n_m + (-1) ** s * math.factorial(n - s) * r ** (n - 2 * s) /   \
                        (math.factorial(s) * math.factorial((n + m) / 2 -s) * math.factorial((n - m) / 2 -s))
                if (m == 0) and (n % 2 == 0):
                    self._polynomials[j_index, :, :] = np.sqrt(n + 1) * r_n_m
                    j_index += 1
                    cos_sin = not(cos_sin)
                else:
                    if cos_sin:
                        self._polynomials[j_index, :, :] = np.sqrt(n + 1) * r_n_m * np.sqrt(2) * np.cos(m * theta)
                        j_index += 1
                        self._polynomials[j_index, :, :] = np.sqrt(n + 1) * r_n_m * np.sqrt(2) * np.sin(m * theta)
                        j_index += 1
                    else:
                        self._polynomials[j_index, :, :] = np.sqrt(n + 1) * r_n_m * np.sqrt(2) * np.sin(m * theta)
                        j_index += 1
                        self._polynomials[j_index, :, :] = np.sqrt(n + 1) * r_n_m * np.sqrt(2) * np.cos(m * theta)
                        j_index += 1

    def __init__(self, *args, **kwargs):

        input_norm = kwargs.pop('normalization', 'unit_rms')
        if input_norm != 'unit_rms':
            raise ValueError('Selected normalization not supported. Noll form only supports unit_rms')

        ZernikeBasis.__init__(self, *args, **kwargs)
        self._normalization = 'unit_rms'

        # Calculate the Noll Zernikes
        self.__calc_noll_zerns__()

        # Noll Zernikes always have unit RMS
        self._rms_factors = np.ones(self._numpolys)

    # def __init__(self, *args, **kwargs):
    #
    #     input_norm = kwargs.pop('normalization', 'unit_rms')
    #     if input_norm not in {'unit_rms'}:
    #         raise ValueError('Only unit_rms is supported by Noll form of Zernikes')
    #     ZernikeBasis.__init__(self, *args, normalization=input_rms, **kwargs)
    #     self._polynomial_map = np.zeros(self.numpolys)
    #
    #     for index in range(0, self.numpolys + 1):
    #         if self._radial_terms[index]:
    #             self._polynomial_map[index] = index - self._term_order[index] / 2
    #         else:
    #             if (self._term_order[index] % 2) == 0:
    #                 # What is the index of the radial term for this order in our standard ordering
    #                 order_start = self._term_order[index] * (self._term_order[index] + 1) / 2
    #                 next_order_start = (self._term_order[index] + 1) * (self._term_order[index] + 2) / 2
    #                 radial_term_index = order_start + self._term_order[index] / 2
    #                 if index < radial_term_index:
    #                     self._polynomial_map[index] = next_order_start - (index - order_start) * 2 - (self._term_order[index] - 1) / 2
    #                 elif index > radial_term_index:
    #                     self._polynomial_map[index] =
    #
    #                 self._polynomial_map[index] =


class FringeZernikesBasisBase(ZernikeBasis):

    def __init__(self, *args, **kwargs):
        ZernikeBasis.__init__(self, *args, **kwargs)

    def show_polynomials(self, plot_terms=None, individual=False):
        if not plot_terms:
            plot_terms = range(0, self._numpolys)

        if individual:
            for polynomial in self._polynomials[plot_terms, :, :]:
                fig, ax = plt.subplots()
                im = ax.imshow(polynomial * self.data_defined, extent=self._extent)
                fig.colorbar(im, ax=ax)
        else:
            min_val = np.min(self._polynomials[:, self._data_defined])
            max_val = np.max(self._polynomials[:, self._data_defined])

            rows = np.floor(np.sqrt(self._numpolys)).astype(int)
            cols = np.ceil(self._numpolys / rows).astype(int)
            plotted = np.zeros((rows, cols), dtype=bool)
            fig, ax = plt.subplots(rows, cols)
            for polynomial, index in zip(self._polynomials[plot_terms, :, :], plot_terms):
                row = np.floor(index / cols).astype(int)
                col = np.mod(index, cols).astype(int)
                im = ax[row, col].imshow(polynomial * self.data_defined, extent=self._extent, vmin=min_val,
                                         vmax=max_val)
                ax[row, col].axis('off')
                plotted[row, col] = True

            for cur_ax in ax[np.nonzero(np.logical_not(plotted))]:
                cur_ax.axis('off')
            fig.subplots_adjust(right=0.825)
            cbar_ax = fig.add_axes([0.85, 0.1, 0.03, 0.8])
            fig.colorbar(im, cax=cbar_ax)


class FringeZernikeBasis(FringeZernikesBasisBase):
    """Class that implements tabulated Fringe Zernikes, as returned by many optical design programs such as Zemax or
        or CodeV."""

    def __init__(self, *args, **kwargs):

        # Set the order needed to get all of the fringe Zernikes
        numpolys = args[0]
        arglist = list(args)
        # TODO:  Put in some way to determine on the max order really needed to avoid unneeded computing here.
        arglist[0] = 12  # Maximum polynomial order required to calculate all the needed Fringe Zernike terms.
        args = tuple(i for i in arglist)
        FringeZernikesBasisBase.__init__(self, *args, **kwargs)

        self._numpolys = numpolys
        self._polynomial_map = [0, 2, 1, 4, 5, 3, 8, 7, 12, 9, 6, 13, 11, 18, 17, 24, 14, 10, 19, 16, 25,
                                23, 32, 31, 40, 20, 15, 26, 22, 33, 30, 41, 39, 50, 49, 60, 84]
        self._polynomials = self._polynomials[self._polynomial_map[0:numpolys], :, :]
        self._term_order = self._term_order[self._polynomial_map[0:numpolys]]
        self._rms_factors = self._rms_factors[self._polynomial_map[0:numpolys]]
        self._radial_terms = self._radial_terms[self._polynomial_map[0:numpolys]]

        # Fringe Zernike are a little weird to state order, since it doesn't have all terms for some of them.
        # However, for the fitting code to work there needs to be a value here, so lets put whatever we're using
        # For the max order
        self._order = arglist[0]


class FringeZernikeBasisHardcoded(FringeZernikesBasisBase):
    """Class that implements tabulated Fringe Zernikes, as returned by many optical design programs such as Zemax or
        or CodeV."""

    def __calc_fringe_zerns__(self):
        self._polynomials = np.zeros((self._numpolys, self._pts_x, self._pts_y,))

        if self._numpolys > 37:
            raise ValueError('More Zernikes requested than defined in Fringe set')

        # Precalculate needed orders of the radial coordinate rho
        # TODO:  Put in some way to determine on the max order really needed to avoid unneeded computing here.
        rho_max_order = 12 # Maximum polynomial order required to calculate all the needed Fringe Zernike terms.
        rho = np.ones((rho_max_order + 1, self._pts_x, self._pts_y))
        rho[1, :, :] = self._r
        for index in range(2, rho_max_order + 1):
            rho[index, :, :] = rho[index - 1, :, :] * rho[1, :, :]

        self._term_order = np.zeros(self._numpolys).astype(int)
        self._radial_terms = np.zeros(self._numpolys).astype(bool)

        # Calculate the polar angle
        theta = np.arctan2(self._ynorm, self._xnorm)

        # Calculate the Zernike value up to self._numpolys.  Since these are tabulated and using if statements to only
        # get the polynomials we want would be a pain, we use a "try" and catch the exception where we have run out of
        # the range of allocated self._polynomials array.
        try:
            self._polynomials[0, :, :] = 1
            self._term_order[0] = 0
            self._radial_terms[0] = True
            self._polynomials[1, :, :] = self._xnorm
            self._term_order[1] = 1
            self._polynomials[2, :, :] = self._ynorm
            self._term_order[2] = 1
            self._polynomials[3, :, :] = 2. * rho[2, :, :] - 1
            self._radial_terms[3] = True
            self._term_order[3] = 2
            self._polynomials[4, :, :] = rho[2, :, :] * np.cos(2 * theta)
            self._term_order[4] = 2
            self._polynomials[5, :, :] = rho[2, :, :] * np.sin(2 * theta)
            self._term_order[5] = 2
            self._polynomials[6, :, :] = (3 * rho[3, :, :] - 2 * rho[1, :, :]) * np.cos(theta)
            self._term_order[6] = 3
            self._polynomials[7, :, :] = (3 * rho[3, :, :] - 2 * rho[1, :, :]) * np.sin(theta)
            self._term_order[7] = 3
            self._polynomials[8, :, :] = 6 * rho[4, :, :] - 6 * rho[2, :, :] + 1
            self._term_order[8] = 4
            self._radial_terms[8] = True
            self._polynomials[9, :, :] = rho[3, :, :] * np.cos(3 * theta)
            self._term_order[9] = 3
            self._polynomials[10, :, :] = rho[3, :, :] * np.sin(3 * theta)
            self._term_order[10] = 3
            self._polynomials[11, :, :] = (4 * rho[4, :, :] - 3 * rho[2, :, :]) * np.cos(2 * theta)
            self._term_order[11] = 4
            self._polynomials[12, :, :] = (4 * rho[4, :, :] - 3 * rho[2, :, :]) * np.sin(2 * theta)
            self._term_order[12] = 4
            self._polynomials[13, :, :] = (10 * rho[5, :, :] - 12 * rho[3, :, :] + 3 * rho[1, :, :]) * np.cos(theta)
            self._term_order[13] = 5
            self._polynomials[14, :, :] = (10 * rho[5, :, :] - 12 * rho[3, :, :] + 3 * rho[1, :, :]) * np.sin(theta)
            self._term_order[14] = 5
            self._polynomials[15, :, :] = 20 * rho[6, :, :] - 30 * rho[4, :, :] + 12 * rho[2, :, :] - 1
            self._term_order[15] = 6
            self._radial_terms[15] = True
            self._polynomials[16, :, :] = rho[4, :, :] * np.cos(4 * theta)
            self._term_order[16] = 4
            self._polynomials[17, :, :] = rho[4, :, :] * np.sin(4 * theta)
            self._term_order[17] = 4
            self._polynomials[18, :, :] = (5 * rho[5, :, :] - 4 * rho[3, :, :]) * np.cos(3 * theta)
            self._term_order[18] = 5
            self._polynomials[19, :, :] = (5 * rho[5, :, :] - 4 * rho[3, :, :]) * np.sin(3 * theta)
            self._term_order[19] = 5
            self._polynomials[20, :, :] = (15 * rho[6, :, :] - 20 * rho[4, :, :] + 6 * rho[2, :, :]) * np.cos(2 * theta)
            self._term_order[20] = 6
            self._polynomials[21, :, :] = (15 * rho[6, :, :] - 20 * rho[4, :, :] + 6 * rho[2, :, :]) * np.sin(2 * theta)
            self._term_order[21] = 6
            self._polynomials[22, :, :] = (35 * rho[7, :, :] - 60 * rho[5, :, :] + 30 * rho[3, :, :] -
                                           4 * rho[1, :, :]) * np.cos(theta)
            self._term_order[22] = 7
            self._polynomials[23, :, :] = (35 * rho[7, :, :] - 60 * rho[5, :, :] + 30 * rho[3, :, :] -
                                           4 * rho[1, :, :]) * np.sin(theta)
            self._term_order[23] = 7
            self._polynomials[24, :, :] = (70 * rho[8, :, :] - 140 * rho[6, :, :] + 90 * rho[4, :, :] -
                                           20 * rho[2, :, :] + 1)
            self._term_order[24] = 8
            self._radial_terms[24] = True
            self._polynomials[25, :, :] = rho[5, :, :] * np.cos(5 * theta)
            self._term_order[25] = 5
            self._polynomials[26, :, :] = rho[5, :, :] * np.sin(5 * theta)
            self._term_order[26] = 5
            self._polynomials[27, :, :] = (6 * rho[6, :, :] - 5 * rho[4, :, :]) * np.cos(4 * theta)
            self._term_order[27] = 6
            self._polynomials[28, :, :] = (6 * rho[6, :, :] - 5 * rho[4, :, :]) * np.sin(4 * theta)
            self._term_order[28] = 6
            self._polynomials[29, :, :] = (21 * rho[7, :, :] - 30 * rho[5, :, :] +
                                           10 * rho[3, :, :]) * np.cos(3 * theta)
            self._term_order[29] = 7
            self._polynomials[30, :, :] = (21 * rho[7, :, :] - 30 * rho[5, :, :] +
                                           10 * rho[3, :, :]) * np.sin(3 * theta)
            self._term_order[30] = 7
            self._polynomials[31, :, :] = (56 * rho[8, :, :] - 105 * rho[6, :, :] + 60 * rho[4, :, :] -
                                           10 * rho[2, :, :]) * np.cos(2 * theta)
            self._term_order[31] = 8
            self._polynomials[32, :, :] = (56 * rho[8, :, :] - 105 * rho[6, :, :] + 60 * rho[4, :, :] -
                                           10 * rho[2, :, :]) * np.sin(2 * theta)
            self._term_order[32] = 8
            self._polynomials[33, :, :] = (126 * rho[9, :, :] - 280 * rho[7, :, :] + 210 * rho[5, :, :] -
                                           60 * rho[3, :, :] + 5 * rho[1, :, :]) * np.cos(theta)
            self._term_order[33] = 9
            self._polynomials[34, :, :] = (126 * rho[9, :, :] - 280 * rho[7, :, :] + 210 * rho[5, :, :] -
                                           60 * rho[3, :, :] + 5 * rho[1, :, :]) * np.sin(theta)
            self._term_order[34] = 9
            self._polynomials[35, :, :] = (252 * rho[10, :, :] - 630 * rho[8, :, :] + 560 * rho[6, :, :] -
                                           210 * rho[4, :, :] + 30 * rho[2, :, :] - 1)
            self._radial_terms[35] = True
            self._term_order[35] = 10
            self._polynomials[36, :, :] = (924 * rho[12, :, :] - 2772 * rho[10, :, :] + 3150 * rho[8, :, :] -
                                           1680 * rho[6, :, :] + 420 * rho[4, :, :] - 42 * rho[2, :, :] + 1)
            self._term_order[36] = 12
            self._radial_terms[36] = True
        except IndexError:
            pass

    def __init__(self, *args, **kwargs):

        FringeZernikesBasisBase.__init__(self, *args, **kwargs)

        self._numpolys = args[0]
        # Fringe Zernike are a little weird to state order, since it doesn't have all terms for some of them.
        self._order = None

        # Calculate the fringe Zernikes
        self.__calc_fringe_zerns__()

        # Figure out the factors needed to give these unit RMS
        self._rms_factors = np.zeros(self._numpolys)
        for index1 in range(0, self._numpolys):
            self._rms_factors[index1] = np.sqrt(np.mean(np.power(self._polynomials[index1][self._data_defined], 2)))


class Legendre2DBasis(PolynomialBasis):
    """Class that implements two-dimensional Legendre polynomial Basis"""

    def __calc_legendres__(self):

        # Calculate the 1D Legendre polynomials of all the orders we're interested in
        poly_x1d = np.zeros((self._order + 1, self._pts_x))
        poly_y1d = np.zeros((self._order + 1, self._pts_y))
        for index in range(0, self._order + 1):
            leg_poly1d = sp.legendre(index)
            poly_x1d[index, :] = leg_poly1d(self._xnorm)
            poly_y1d[index, :] = leg_poly1d(self._ynorm).T

        # Calculate the 2D polynomials by taking the product of each of the x, y 1D polynomials
        # self._numpolys = (self._order + 1) ** 2
        self._numpolys = int(((self._order + 1) * (self._order + 2)) / 2)
        polys_2d = np.zeros((self._order + 1, self._order + 1, self._pts_x, self._pts_y))
        for index_i in range(0, self._order + 1):
            for index_j in range(0, self._order + 1):
                polys_2d[index_i, index_j, :, :] = np.einsum('i, j->ij', poly_x1d[index_i, :], poly_y1d[index_j, :])
        count = 0

        # Reorder the polynomials according to order.  We won't use the ones that don't have all the terms.  In other
        # words make a Zernike like pyramid of terms, but with a single index ordering.
        self._polynomials = np.zeros((self._numpolys, self._pts_x, self._pts_y))
        self._term_order = np.zeros(self._numpolys, dtype=int)
        for index_i in range(0, self._order + 1):
            for index_j in range(index_i, -1, -1):
                self._term_order[count] = index_i
                self._polynomials[count, :, :] = polys_2d[index_j, index_i - index_j, :, :]
                count += 1

    def __init__(self, *args, **kwargs):

        input_norm = kwargs.pop('normalization', 'natural')
        if input_norm not in {'natural', 'unit_pv', 'unit_rms'}:
            raise ValueError('Selected normalization not supported')
        PolynomialBasis.__init__(self, *args, normalization=input_norm, **kwargs)

        # Calculate the legendres
        self.__calc_legendres__()

        self._data_defined = (np.abs(self._xnorm) <= 1) & (np.abs(self._ynorm) <= 1)

        # Figure out the factors needed to give these unit RMS
        self._rms_factors = np.zeros(self._numpolys)
        for index1 in range(0, self._numpolys):
            self._rms_factors[index1] = np.sqrt(np.mean(np.power(self._polynomials[index1][np.transpose(self._data_defined)], 2)))
            # plt.figure()
            # plt.imshow(self._polynomials[index1])
            # plt.show()

class GramSchmidtBasis(PolynomialBasis):
    """Class that implements a Gram-Schmidt orthonormalized version of the the polynomials given in the second argument.
        Output polynomial should be exactly orthogonal over the boolean array defining the points in the array where the
        data is defined."""

    def __init__(self, *args, **kwargs):

        new_data_defined = args[0]
        nonorth_polynomials = args[1]

        input_norm = kwargs.pop('normalization', 'natural')
        if input_norm not in {'natural', 'unit_pv', 'unit_rms'}:
            raise ValueError('Selected normalization not supported')

        PolynomialBasis.__init__(self, nonorth_polynomials.order, nonorth_polynomials.pts_x, nonorth_polynomials.pts_y,
                                 normalization=input_norm, **kwargs)

        self.numpolys = nonorth_polynomials.numpolys
        self._term_order = nonorth_polynomials._term_order
        self._data_defined = new_data_defined

        self._polynomials = np.zeros((self._numpolys, self._xnorm.shape[1], self._ynorm.shape[0]))
        self._rms_factors = np.zeros(self.numpolys)
        self._rms_factors[0] = 1

        #  Do orthonormalization
        # First polynomial always stay the same
        self._polynomials[0, :, :] = nonorth_polynomials.polynomials[0, :, :]
        # Loop through all of the polynomials that need to be orthogonalized, keeping track of an index
        # for accessing other polynomials as well
        self._d = np.zeros((nonorth_polynomials.numpolys, self.numpolys))
        for init_polynomial, index1 in zip(nonorth_polynomials.polynomials[1:, :, :],
                                           range(1, nonorth_polynomials.numpolys)):
            # Start out new polynomial by setting it to the current old polynomial's value
            self._polynomials[index1, :, :] = np.copy(init_polynomial)
            # Loop over all defined polynomials up to current order, using their contributions to orthogonalize
            # current order
            for intermed_polynomial, index2 in zip(self._polynomials[0:index1, :, :],
                                                   range(0, index1)):
                # Calculate the projection of lower order polynomials on current one
                cur_d = -np.einsum('i, i', init_polynomial[new_data_defined], intermed_polynomial[new_data_defined]) \
                        / np.einsum('i, i', intermed_polynomial[new_data_defined],
                                    intermed_polynomial[new_data_defined])
                # Subtract the lower order polynomial's contribution to current one
                self._polynomials[index1, :, :] += cur_d * intermed_polynomial
                # Store the coefficients from this projection for use in fitting code.
                self._d[index1, index2] = cur_d

            self._rms_factors[index1] = np.sqrt(np.mean(np.power(self._polynomials[index1][new_data_defined], 2)))

            # Mask out region over which new polynomials are orthogonal
            self._polynomials[index1, :, :] = self._polynomials[index1, :, :] * new_data_defined


# @property
# def polynomials(self):
#     if self._normalization == 'natural':
#         return self._polynomials
#     elif self._normalization == 'unit_pv':
#         scale = [1]
#         offset = [0]
#         for polynomial in self._polynomials[1:, :, :]
#             max_val = np.max(polynomial[self.data_defined])
#             min_val = np.min(polynomial[self.data_defined])
#             scale.append(2 / (max_val - min_val))
#             offset.append(-(min_val * scale + 1))
#         return (np.einsum('ijk, i->ijk', self._polynomials, scale)) + offset
#     elif self._normalization == 'unit_rms':
#         return (np.einsum('ijk, i->ijk', self._polynomials, 1 / self._rms_factors))


class Wavefront:
    """Base class for any wavefront.  Wavefronts can be represented by polynomials and coefficients or can be
    coefficients of a basis set, for which an instance of the class needs to be included"""

    def __init__(self, *args, **kwargs):
        self._data_diameter = args[0]
        self._wavelength = args[1]
        self._wf_units = args[2]

        self._data = None
        self._data_defined = None
        self._basis = None

        # # Calculate positive and negative maximum number of points to have
        # min_x = -np.fix(self._pts_x / 2)
        # max_x = np.fix((self._pts_x + 1) / 2)
        # min_y = -np.fix(self._pts_y / 2)
        # max_y = np.fix((self._pts_y + 1) / 2)
        #
        # # Step sizes to get -1 to 1 range
        # self._step_x = 1 / np.fix((self._pts_x - 1) / 2)
        # self._step_y = 1 / np.fix((self._pts_y - 1) / 2)
        #
        # # Create open arrays for coordinates and scales
        # self._ynorm, self._xnorm = np.ogrid[min_y:max_y, min_x:max_x]
        # self._xnorm = self._xnorm * self._step_x
        # self._ynorm = self._ynorm * self._step_y
        #
        # self._extent = [min_x * self._step_x, max_x * self._step_x, min_y * self._step_y, max_y * self._step_y]

        """Base class constructor that sets attributes that every wavefront will need.  At the time of writing the
            was that there would be no common attributes between a sampled, point-by-point, representation and a basis
            polynomial representation.  As a consequence this constructor is a placeholder and there are no attributes
       Arguments:
            None"""

    @property
    def pv(self):
        minval = np.min(self.data * self.data_defined)
        maxval = np.max(self.data * self.data_defined)
        return (maxval - minval) * self._wf_units

    @property
    def rms(self):
        rms = np.sqrt(np.mean(np.power(self._data[self.data_defined], 2)))
        return rms * self._wf_units

    @property
    def data_diameter(self):
        return self._data_diameter

    @data_diameter.setter
    def data_diameter(self, new_data_diameter):
        self._data_diameter = new_data_diameter


class PolynomialWavefront(Wavefront):
    """Class that represents a wavefront represented by a set of coefficients and an instance of the class containing
        the polynomial set used, a descendent of the PolynomialBasis class"""

    def __init__(self, *args, **kwargs):

        Wavefront.__init__(self, *args[2:], **kwargs)
        self._coeffs = args[0]  # List of wavefront coefficients, has units of length or phase.
        # self._data_diameter = args[2] #The physical width of the wavefront represented by this array

        self._basis = args[1]  # Polynomial basis with numpolys = number of elements of wavefront_coeffs

    @property
    def data(self):
        return np.einsum('i,ijk->jk', self._coeffs, self._basis.polynomials)

    @property
    def data_defined(self):
        return self._basis.data_defined

    @data_defined.setter
    def data_defined(self, new_data_defined, new_basis_instance=False):
        if new_basis_instance:
            self._basis = copy.deepcopy(self._basis)
        self._basis._data_defined = new_data_defined

    @property
    def coeffs(self):
        return self._coeffs

    @coeffs.setter
    def coeffs(self, new_coeffs):
        self._coeffs = new_coeffs

    @property
    def rms(self):
        return np.sqrt(np.einsum('i->', (self._coeffs * self._basis.rms_factors) ** 2)) * self._wf_units

    @property
    def pts_x(self):
        return self._basis.pts_x

    @property
    def pts_y(self):
        return self._basis.pts_y

    @property
    def pts(self):
        return self._basis.pts_x, self._basis.pts_y

    @pts_x.setter
    def pts_x(self, new_pts_x, new_basis_instance=False):
        if new_basis_instance:
            self._basis = copy.deepcopy(self._basis)
        self._basis.pts_x = new_pts_x

    @pts_y.setter
    def pts_y(self, new_pts_y, new_basis_instance=False):
        if new_basis_instance:
            self._basis = copy.deepcopy(self._basis)
        self._basis.pts_y = new_pts_y

    @pts.setter
    def pts(self, new_pts, new_basis_instance=False):
        if new_basis_instance:
            self._basis = copy.deepcopy(self._basis)
        self._basis.pts_x = new_pts[0]
        self._basis.pts_y = new_pts[1]

    @property
    def order(self):
        return self._basis.order

    @order.setter
    def order(self, new_order, new_basis_instance=False):
        if new_basis_instance:
            self._basis = copy.deepcopy(self._basis)
        self._basis.order = new_order

    @property
    def step_x(self):
        return self._basis.step_x * self._data_diameter / 2

    @property
    def step_y(self):
        return self._basis.step_y * self._data_diameter / 2

    @property
    def extent(self):
        return self._basis.extent * self._data_diameter / 2

    def show_wavefront(self):
        fig, ax = plt.subplots()
        curextent = self.extent.to(units.mm).magnitude
        im = ax.imshow(self.data * self.data_defined, extent=curextent)
        fig.colorbar(im, ax=ax)


class PointByPointWavefront(Wavefront):
    """Class that represents a wavefront on a regularly spaced grid of values at individual points"""

    @property
    def pts_x(self):
        return self._pts_x

    @property
    def pts_y(self):
        return self._pts_y

    @property
    def pts(self):
        return self._pts_x, self._pts_y

    @pts_x.setter
    def pts_x(self, new_pts_x, scaling='scale'):
        if scaling == 'fixed':
            if self._pts_x > new_pts_x:
                em_vec = embed(self._pts_x, new_pts_x)
                self._data = self._data[em_vec, :]
            else:
                em_vec = embed(new_pts_x, self._pts_x)
                data_new = np.zeros((new_pts_x, self._pts_y))
                data_new[em_vec, :] = self._data
                self._data = data_new
        elif scaling == 'scale':
            print("not implemented")
        self._pts_x = new_pts_x

    @pts_y.setter
    def pts_y(self, new_pts_y, ):
        self._pts_y = new_pts_y

    @pts.setter
    def pts(self, new_pts):
        self._pts_x = new_pts[0]
        self._pts_y = new_pts[1]

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, new_data):
        self._data = new_data

    @property
    def data_defined(self):
        return self._data_defined

    @data_defined.setter
    def data_defined(self, new_data_defined):
        self._data_defined = new_data_defined

    @property
    def extent(self):
        radius = self.data_diameter / 2
        radius = radius.to(units.mm).magnitude
        cur_extent = np.array([-radius, radius, -radius, radius]) * units.mm
        return cur_extent
        # return np.array([-self.data_diameter/2, self.data_diameter/2, -self.data_diameter/2, self.data_diameter/2])

    def __init__(self, *args, **kwargs):
        Wavefront.__init__(self, *args[2:], **kwargs)
        self._data = args[0]  # Values of phase/opd in a 2D array
        pts = self._data.shape
        self._pts_x = pts[0]
        self._pts_y = pts[1]
        self._data_defined = args[1]  # 2D array of booleans, True where valid data exists.

        if 'wf_basis' in kwargs:
            self._basis = kwargs['wf_basis']
        else:
            self._basis = None

    @property
    def coeffs(self, wf_basis=None):
        if self._basis is None:
            if wf_basis is None:
                raise ValueError('Must provide a basis first to get coefficients')
            else:
                self._basis = wf_basis
        elif wf_basis is not None:
            self._basis = wf_basis

        return self._basis.polyfit(self._data, self._data_defined)

    def show_wavefront(self):
        fig, ax = plt.subplots()
        curextent = self.extent.to(units.mm).magnitude
        im = ax.imshow(self._data * self.data_defined, extent=curextent)
        fig.colorbar(im, ax=ax)


def main():
    # a = ZernikeBasis(4, 16, 16)
    a = NollZernikeBasis(7, 1024, 1024)
    a.show_polynomials()

    a = Legendre2DBasis(6, 128, 128)
    a.show_polynomials()

    # new_data_defined = (abs(a.xnorm) <= 0.7) & (abs(a.ynorm) <= 0.7)
    # new_data_defined = np.sqrt(a.xnorm ** 2 + a.ynorm ** 2) <= 0.7

    # b = GramSchmidtBasis(new_data_defined, a)

    # coeffs = np.zeros(a.numpolys)
    # coeffs[1] = 2
    coeffs = np.random.randn(a.numpolys)
    # coeffs = np.ones(a.numpolys)
    print(coeffs)

    c = PolynomialWavefront(coeffs, a, 5 * units.mm, 633, units.nm)

    # b = Legendre2DBasis(int(np.max(a._term_order)), 128, 128)
    # b.show_polynomials()

    # q = FringeZernikeBasisHardcoded(5, 128, 128)
    # q.show_polynomials()

    # x = FringeZernikeBasis(5, 128, 128)
    # x.show_polynomials()

    # TODO: Clean up stuff involving units and Matplotlib/imshow and extent
    #       Think through how to set extent and generally track spatial scale
    #       Complete scaling/interpolation code for point-by-point
    #       Also, implement some version of Hexikes.  Could be just tabulated for now
    #       Add properties/setters for Legendre and Fringe zernikes.
    #       Try to make sure different basis sets are consistent with x/y flipping, etc.
    #       Decide whether to leave hardcoded Fringe Zernikes in.  Decide what to do with
    #       Plotting function for Fringe Zernikes.

    # d = PointByPointWavefront(c.data, c.data_defined, 3 * units.mm, 633, units.nm, wf_basis=b)

    # e = PolynomialWavefront(d.coeffs, b, 3 * units.mm, 633, units.nm)

    # print(d.coeffs)

    print(c.pv)
    print(c.rms)
    # print(d.pv)
    # print(d.rms)

    # e.data_defined = d.data_defined
    c.show_wavefront()
    # d.show_wavefront()
    # e.show_wavefront()
    # fig, ax = plt.subplots()
    # im = ax.imshow(c.data * a.data_defined, extent=a.extent)
    # fig.colorbar(im, ax=ax)

    # fig, ax = plt.subplots()
    # im = ax.imshow(d.data * a.data_defined, extent=a.extent)
    # fig.colorbar(im, ax=ax)
    #
    # fig, ax = plt.subplots()
    # im = ax.imshow(e.data * a.data_defined, extent=a.extent)
    # fig.colorbar(im, ax=ax)

    plt.show()

    # print(b.rms_factors)
    #
    # for zernike in b.polynomials:
    #     rms = np.sqrt(np.mean(np.power(zernike[b.data_defined], 2)))
    #     # cross = np.einsum('i, i', b.polynomials[7][new_data_defined], zernike[new_data_defined]) / \
    #     #        np.einsum('i, i', zernike[new_data_defined], zernike[new_data_defined])
    #     print(rms)
    #     fig, ax = plt.subplots()
    #     im = ax.imshow(zernike * b.data_defined, extent=b.extent)
    #     fig.colorbar(im, ax=ax)
    #     plt.show()

    # for zernike in a.polynomials:
    #     rms = np.sqrt(np.mean(np.power(zernike[a.data_defined], 2)))
    #     # cross = np.einsum('i, i', b.polynomials[7][new_data_defined], zernike[new_data_defined]) / \
    #     #        np.einsum('i, i', zernike[new_data_defined], zernike[new_data_defined])
    #     print(rms)
    #     fig, ax = plt.subplots()
    #     im = ax.imshow(zernike * a.data_defined, extent=a.extent)
    #     fig.colorbar(im, ax=ax)
    #     plt.show()


if __name__ == '__main__':
    # units = pint.UnitRegistry()
    units.setup_matplotlib(True)
    # units.load_definitions('pint_waves_def.txt')
    units.define('waves = radian * 2 * pi =  wave = wv')
    # units.define('waves = [] =  wave = wv')
    wavelength = 632.8 * units.nm
    c = pint.Context('optics')
    c.add_transformation('meter', 'radian',
                         lambda units, x: x * 2 * np.pi / wavelength)
    c.add_transformation('radian', 'meter',
                         lambda units, x: x * wavelength / (2 * np.pi))
    units.add_context(c)
    units.enable_contexts('optics')

    myval = 1 * units.wv
    # print(myval.to(units.nm))
    # print(myval.to(units.rad))

    units.disable_contexts()
    units.remove_context('optics')
    del (c)
    c = pint.Context('optics')
    # c.remove_transformation('meter', 'radian')
    # c.remove_transformation('radian', 'meter')
    wavelength = 500 * units.nm
    c.add_transformation('meter', 'radian',
                         lambda units, x: x * 2 * np.pi / wavelength)
    c.add_transformation('radian', 'meter',
                         lambda units, x: x * wavelength / (2 * np.pi))
    units.add_context(c)
    units.enable_contexts('optics')
    # print(myval.to(units.nm))
    # print(myval.to(units.rad))

    main()
