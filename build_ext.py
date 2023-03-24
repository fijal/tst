
from pathlib import Path
from cffi import FFI

ffibuilder = FFI()
header = open(Path(__file__).parent / "sim_tb.h", "r").read()
header.split("\n")
ffibuilder.cdef(header + "\n" + """
void* malloc(size_t);
""")
ffibuilder.set_source("_sim_c", '#include "sim_tb.h"', sources=["c_source.c"])

if __name__ == '__main__':
    ffibuilder.compile(verbose=True)
