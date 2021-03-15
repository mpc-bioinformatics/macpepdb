import signal

from .. import process_context

class GenericProcess(process_context.Process):
    """
    Generic process, inherit from the applications process context. Ignores CTRT-C/SIGINT
    """
    def __init__(self):
        super().__init__()
        # Ignore CTRL-C/SIGINT
        signal.signal(signal.SIGINT, signal.SIG_IGN)
