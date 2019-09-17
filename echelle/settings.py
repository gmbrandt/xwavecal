import operator

from banzai.utils.file_utils import ccdsum_to_filename
from banzai.settings import make_calibration_filename_function
from banzai.utils.instrument_utils import InstrumentCriterion

from banzai import settings

from banzai_nres.utils.fiber_utils import fibers_state_to_filename
from banzai_nres.utils.runtime_utils import get_telescope_filename


settings.FRAME_CLASS = 'banzai_nres.images.NRESImage'

settings.FRAME_SELECTION_CRITERIA = [InstrumentCriterion('type', operator.contains, 'NRES')]

settings.ORDERED_STAGES = ['banzai.bpm.BPMUpdater',
                           'banzai.qc.SaturationTest',
                           'banzai.bias.OverscanSubtractor',
                           'banzai.gain.GainNormalizer',
                           'banzai.trim.Trimmer',
                           'banzai.bias.BiasSubtractor',
                           'banzai.dark.DarkSubtractor',
                           'banzai_nres.traces.LoadTrace',
                           'banzai_nres.extract.RectifyTwodSpectrum',
                           'banzai_nres.extract.BoxExtract']

settings.CALIBRATION_MIN_FRAMES = {'BIAS': 5,
                                   'DARK': 3,
                                   'LAMPFLAT': 5,
                                   'TRACE': 1}

# Trace settings #
TRACE_FIT_INITIAL_DEGREE_TWO_GUESS = 90  # DO NOT HAPHAZARDLY CHANGE THIS
TRACE_FIT_POLYNOMIAL_ORDER = 4  # DO NOT HAPHAZARDLY CHANGE THIS
TRACE_TABLE_NAME = 'TRACE'
WINDOW_FOR_TRACE_IDENTIFICATION = {'max': 2100, 'min': 2000}  # pixels
MIN_FIBER_TO_FIBER_SPACING = 10  # pixels
MIN_SNR_FOR_TRACE_IDENTIFICATION = 6
# Blaze settings.
BLAZE_TABLE_NAME = 'BLAZE'
BLAZE_CORRECTED_BOX_SPECTRUM_NAME = 'BLZCBOX'
# Extraction settings #
MAX_EXTRACTION_HALF_WINDOW = 10  # pixels
BOX_EXTRACTION_HALF_WINDOW = 10  # pixels
BOX_SPECTRUM_NAME = 'SPECBOX'

# order identification settings #
ref_id = 26  # reference id of the central diffraction order in the database arc template.

# Wavelength calibration settings #
OVERLAP_TABLE_NAME = 'OVERLAP'
EMISSION_LINES_TABLE_NAME = 'LINES'
# Note that we can also accomplish the above by flipping the image.
max_red_overlap = 1000
max_blue_overlap = 2000
global_scale_range = (0.5, 1.5)
overlap_linear_scale_range = (0.5, 2)
approx_detector_range_angstroms = 5000
approx_num_orders = 67
# principle order number (m0) settings for the IdentifyPrincipleOrderNumber stage
m0_range = (50, 54)  # start (inclusive), stop (exclusive)
principle_order_number = 52
min_num_overlaps = 5 # the minimum number of overlaps which must be fit for the wavelength solution to proceed.
# refine wavelength solution settings
min_peak_snr = 10  # the min signal to noise for an emission peak to be considered in the wavelength solution.
# initial model that is constrained via the overlaps and is used to find the global scale:
initial_wavelength_model = {1: [0, 1, 2],
                            2: [0, 1, 2]}
# wavelength model for the initial refine stage:
intermediate_wavelength_model = {0: [0, 1, 2],
                                 1: [0, 1, 2],
                                 2: [0, 1, 2]}
# wavelength model that the final refine stage will end at (stage begins with the intermediate model):
final_wavelength_model = {0: [0, 1, 2, 3, 4, 5],
                          1: [0, 1, 2, 3, 4, 5],
                          2: [0, 1, 2, 3, 4, 5],
                          3: [0, 1, 2, 3, 4, 5],
                          4: [0]}

# general reduction settings:

settings.CALIBRATION_SET_CRITERIA = {'BIAS': ['ccdsum'],
                                     'DARK': ['ccdsum'],
                                     'LAMPFLAT': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit'],
                                     'TRACE': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit'],
                                     'BLAZE': ['ccdsum', 'fiber0_lit', 'fiber1_lit', 'fiber2_lit']}

settings.CALIBRATION_FILENAME_FUNCTIONS = {'BIAS': make_calibration_filename_function('BIAS', [ccdsum_to_filename],
                                                                                      get_telescope_filename),
                                           'DARK': make_calibration_filename_function('DARK', [ccdsum_to_filename],
                                                                                      get_telescope_filename),
                                           'LAMPFLAT': make_calibration_filename_function('LAMPFLAT',
                                                                                          [ccdsum_to_filename,
                                                                                           fibers_state_to_filename],
                                                                                          get_telescope_filename),
                                           'TRACE': make_calibration_filename_function('TRACE', [ccdsum_to_filename,
                                                                                       fibers_state_to_filename],
                                                                                       get_telescope_filename)}

settings.CALIBRATION_IMAGE_TYPES = ['BIAS', 'DARK', 'LAMPFLAT', 'TRACE']

settings.LAST_STAGE = {'BIAS': 'banzai.trim.Trimmer',
                       'DARK': 'banzai.bias.BiasSubtractor',
                       'LAMPFLAT': 'banzai.dark.DarkSubtractor',
                       'TRACE': 'banzai.dark.DarkSubtractor',
                       'DOUBLE': None,
                       'TARGET': None}

settings.EXTRA_STAGES = {'BIAS': ['banzai.bias.BiasMasterLevelSubtractor', 'banzai.bias.BiasComparer'],
                         'DARK': ['banzai.dark.DarkNormalizer', 'banzai.dark.DarkComparer'],
                         'LAMPFLAT': None,
                         'TRACE': None,
                         'DOUBLE': ['banzai_nres.fibers.IdentifyFibers',
                                    'banzai_nres.blaze.BackgroundSubtractSpectrum',
                                    'banzai_nres.blaze.ApplyBlaze',
                                    'banzai_nres.wavelength.Initialize',
                                    'banzai_nres.wavelength.AddWavelengthColumn',
                                    'banzai_nres.wavelength.LoadReferenceLineList',
                                    'banzai_nres.wavelength.IdentifyArcEmissionLinesLowSN',
                                    'banzai_nres.wavelength.BlazeCorrectArcEmissionLines',
                                    'banzai_nres.wavelength.FitOverlaps',
                                    'banzai_nres.wavelength.IdentifyArcEmissionLines',
                                    #'banzai_nres.wavelength.IdentifyPrincipleOrderNumber',
                                    'banzai_nres.wavelength.SolveFromOverlaps',
                                    'banzai_nres.wavelength.FindGlobalScale',
                                    'banzai_nres.wavelength.SolutionRefineInitial',
                                    'banzai_nres.wavelength.SolutionRefineFinal',
                                    'banzai_nres.wavelength.IdentifyArcEmissionLinesLowSN',
                                    'banzai_nres.wavelength.ApplyToSpectrum',
                                    #'banzai_nres.utils.temp_overlaps.FindOverlapsFromWavelengths',
                                    'banzai_nres.wavelength.TabulateArcEmissionLines'],
                         'TARGET': None}

settings.CALIBRATION_STACKER_STAGE = {'BIAS': 'banzai.bias.BiasMaker',
                                      'DARK': 'banzai.dark.DarkMaker',
                                      'LAMPFLAT': 'banzai_nres.flats.FlatStacker',
                                      'TRACE': 'banzai_nres.traces.TraceMaker'}
