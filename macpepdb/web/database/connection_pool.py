"""
Implements connection pooling for the database.
"""
# std imports
from threading import Semaphore

# external imports
from psycopg2.pool import ThreadedConnectionPool

class ConnectionPool(ThreadedConnectionPool):
    """
    Extends psycopg2.pool.ThreadedConnectionPool by a semaphore,
    forcing a thread to wait if no connection is available.

    Parameters
    ----------
    minconn : int
        Minimum amount of connections
    maxconn : int
        Maximum amount of connections
    """

    def __init__(self, minconn, maxconn, *args, **kwargs):
        self._semaphore = Semaphore(maxconn)
        super().__init__(minconn, maxconn, *args, **kwargs)

    def getconn(self, *args, **kwargs):
        self._semaphore.acquire()
        return super().getconn(*args, **kwargs)

    def putconn(self, *args, **kwargs):
        super().putconn(*args, **kwargs)
        self._semaphore.release()