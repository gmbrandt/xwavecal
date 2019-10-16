import os
import tempfile
import shutil


class Translator(dict):
    def __init__(self, header_keys, type_keys):
        self.header_keys = header_keys
        self.type_keys = type_keys

    def __getitem__(self, key):
        return self.header_keys.get(key, key)


def writeto(hdu_list, filepath, fpack=False, overwrite=True, output_verify='fix+warn', quant=int(1E6)):
    if not fpack:
        hdu_list.writeto(filepath.split('.fz')[0], overwrite=overwrite, output_verify=output_verify)
    else:
        if not filepath.endswith('.fz'):
            filepath += '.fz'
        base_filename = os.path.basename(filepath).split('.fz')[0]
        with tempfile.TemporaryDirectory() as temp_directory:
            hdu_list.writeto(os.path.join(temp_directory, base_filename),
                             overwrite=overwrite, output_verify=output_verify)
            command = 'fpack -q {quantization} {temp_directory}/{basename}'
            os.system(command.format(quantization=int(quant), temp_directory=temp_directory, basename=base_filename))
            base_filename += '.fz'
            shutil.move(os.path.join(temp_directory, base_filename), filepath)
