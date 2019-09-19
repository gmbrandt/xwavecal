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
min_num_overlaps = 5  # the minimum number of overlaps which must be fit for the wavelength solution to proceed.
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

