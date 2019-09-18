import abc


class Stage(object):
    def __init__(self, runtime_context):
        self.runtime_context = runtime_context

    @abc.abstractmethod
    def do_stage(self, image):
        return image
