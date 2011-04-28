
pkg = 'webbpsf'

setupargs = {
    'version'       :      	"0.2",
    'description'   :       'James Webb Space Telescope synthetic PSFs',
    'fullname'      :       'AstroLib WebbPSF',
    'license'       :       'BSD',
    'author'        :     	"Marshall Perrin",
    'author_email'  :      	"mperrin@stsci.edu",
    'url'           :  		"http://www.stsci.edu/~mperrin/software/webbpsf",
    'platforms'     :      	["Linux","Solaris","Mac OS X", "Win"],
    'requires'      :       ['pyfits','numpy', 'matplotlib', 'scipy', 'atpy'],
    'data_files'    :     	[ ( pkg+'/data', [ 'data/generic/*', 'data/wavecat/*' ] ),
                                ( pkg+'/tests', [ 'tests/*'] ),
                            ]
    }

