

class Error(Exception):
    """Base class for exceptions in this module."""
    pass

class ConvergenceError(Error):
    #Exception raised when lammps do not converge.


    def __init__(self, message):
        self.message = message