import matplotlib.pyplot as plt
from webbpsf import roman

#### Create webbpsf-roman_page_header.png
wfi = roman.WFI()

# UNRESOLVED: GRISM0 errors out. poppy's Instrument._get_weights()
# drops wavelengths with throughputs <0.4. GRISM0's peak is well
# below 0.1 and numpy won't take the min/max of an empty array.
filters_no_grism0 = [f for f in wfi.filter_list if f != 'GRISM0']

long = 5 # should be 6 if GRISM0 is included
wide = 2

fig, axs = plt.subplots(wide, long, figsize=(12, 6), sharey=True)

for i, filter in enumerate(sorted(filters_no_grism0)):
    r = int(np.floor(i / long))
    c = i % long
    ax = axs[r][c]

    wfi.filter = filter
    psf = wfi.calc_psf(oversample=4)

    display_psf(psf, ax=ax, colorbar=False, title=filter)
    ax.title.set_fontsize(20)
    ax.tick_params(axis='both', labelsize=10)
    ax.xaxis.label.set_visible(False)
    ax.yaxis.label.set_visible(False)


#axs[-1][-1].remove() # uncomment if GRISM0 is included again

fig.tight_layout(w_pad=.1, h_pad=0)
fig.tight_layout(w_pad=.1, h_pad=0) # calling twice somehow tightens h_pad
#fig.savefig('webbpsf-roman_page_header.png', dpi=100, facecolor='w')

#### Create compare_wfi_sca09_sca17.png

wfi2 = roman.WFI()
wfi2.filter = 'F129'
wfi2.detector = 'SCA09'
wfi2.detector_position = (4, 4)
psf_sca09 = wfi2.calc_psf()
wfi2.detector = 'SCA17'
wfi2.detector_position = (4092, 4092)
psf_sca17 = wfi2.calc_psf()

fig2, (ax_sca09, ax_sca17, ax_diff) = plt.subplots(1, 3, figsize=(16, 4))

webbpsf.display_psf(psf_sca09, ax=ax_sca09, imagecrop=2.0,
                    title='WFI SCA09, bottom left - F129')
webbpsf.display_psf(psf_sca17, ax=ax_sca17, imagecrop=2.0,
                    title='WFI SCA17, top right - F129')
webbpsf.display_psf_difference(psf_sca09, psf_sca17, ax=ax_diff,
                               vmax=5e-3, title='SCA09 - SCA17', imagecrop=2.0)
fig2.tight_layout(w_pad=.5)
#fig2.savefig('compare_wfi_sca09_sca17.png', dpi=100, facecolor='w')
