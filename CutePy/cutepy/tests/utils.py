class DummyComm(object):
    def __init__(self):
        pass

    def Get_rank(self):
        return 0

    def Get_size(self):
        return 1

    def send(self, dest=0):
        pass

    def recv(self, source=0):
        return 0.0


def get_comm_world():
    try:
        from mpi4py import MPI
        comm = MPI.COMM_WORLD
    except ModuleNotFoundError:
        comm = DummyComm()
    return comm
