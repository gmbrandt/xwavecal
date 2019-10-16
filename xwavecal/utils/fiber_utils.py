import numpy as np


def fiber_states_from_header(fiber_string):
    parsed_objects = fiber_string.split('&')
    fiber0_lit = parsed_objects[0].lower() != 'none'
    fiber1_lit = parsed_objects[1].lower() != 'none'
    fiber2_lit = parsed_objects[2].lower() != 'none'
    return int(fiber0_lit), int(fiber1_lit), int(fiber2_lit)


def wavecal_fibers_from_header(fiber_string, wavecal_key='thar'):
    parsed_objects = fiber_string.split('&')
    fiber0_wavecal = parsed_objects[0].lower() == wavecal_key
    fiber1_wavecal = parsed_objects[1].lower() == wavecal_key
    fiber2_wavecal = parsed_objects[2].lower() == wavecal_key
    return int(fiber0_wavecal), int(fiber1_wavecal), int(fiber2_wavecal)


def fibers_state_to_filename(image):
    return str(int(image.fiber0_lit)) + str(int(image.fiber1_lit)) + str(int(image.fiber2_lit))


def lit_wavecal_fibers(image):
    return np.array([0] * image.fiber0_wavecal + [1] * image.fiber1_wavecal + [2] * image.fiber2_wavecal)


def lit_fibers(image):
    return np.array([0] * image.fiber0_lit + [1] * image.fiber1_lit + [2] * image.fiber2_lit)
