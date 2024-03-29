# std imports
import signal

# internal imports
from macpepdb import process_context

class GenericProcess(process_context.Process):
    """
    Generic process, inherit from the applications process context.
    """
    def __init__(self, termination_event: process_context.Event):
        super().__init__()
        self.__termination_event = termination_event
        
    @property
    def termination_event(self) -> process_context.Event:
        """
        Returns
        -------
        Event for calling termination.
        """
        return self.__termination_event

    def activate_signal_handling(self):
        """
        Handles SIGINT and SIGTERM
        """
        signal.signal(signal.SIGINT, self.stop_signal_handler)
        signal.signal(signal.SIGTERM, self.stop_signal_handler)

    def stop_signal_handler(self, signal_number, frame):
        """
        Sets the termination event. Callback for `activate_signal_handling`

        Argument
        ========
        signal_number
            Signal number
        frame
            Signal frame
        """
        self.termination_event.set()