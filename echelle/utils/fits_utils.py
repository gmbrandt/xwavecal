import os
import tempfile
import shutil

import logging




def writeto(hdu_list, filepath, fpack=False, overwrite=True, output_verify='fix+warn'):
    if not fpack:
        hdu_list.writeto(filepath.split('.fz')[0],
                         overwrite=overwrite, output_verify=output_verify)
    else:
        if not filepath.endswith('.fz'):
            filepath += '.fz'
        logger.info('Writing file to {filepath}'.format(filepath=filepath))
        base_filename = os.path.basename(filepath).split('.fz')[0]
        with tempfile.TemporaryDirectory() as temp_directory:
            hdu_list.writeto(os.path.join(temp_directory, base_filename),
                             overwrite=overwrite, output_verify=output_verify)
            if os.path.exists(filepath):
                os.remove(filepath)
            command = 'fpack -q 64 {temp_directory}/{basename}'
            os.system(command.format(temp_directory=temp_directory, basename=base_filename))
            base_filename += '.fz'
            shutil.move(os.path.join(temp_directory, base_filename), filepath)
