import numpy as np

import xwavecal.munge as mg


class TestRotations:
    def test_rotate(self):
        a = np.array([[1, 2], [4, 3]])
        image = type('Image', (), {'data': a, 'ivar': a})
        assert np.allclose(mg.Rot90(None).do_stage(image).data, [[2, 3], [1, 4]])

    def test_flip(self):
        a = (np.arange(2) * np.ones((2, 2))).T
        image = type('Image', (), {'data': a, 'ivar': a})
        assert np.allclose(mg.FlipVert(None).do_stage(image).data, [[1, 1], [0, 0]])
        image.data = np.arange(2) * np.ones((2, 2))
        assert np.allclose(mg.FlipHoriz(None).do_stage(image).data, [[1, 0], [1, 0]])
