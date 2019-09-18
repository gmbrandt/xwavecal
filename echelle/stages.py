import abc


class Stage(object):
    def __init__(self):
        pass

    @abc.abstractmethod
    def do_stage(self, data):
        return data
