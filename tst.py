
import struct, time
import numpy as np
from _sim_c import lib, ffi

SIZE = 16482

from betse.science.chemistry.networks import Channel, MasterOfNetworks

class Cells(object):
    def __init__(self):
        self.mem_i = np.arange(SIZE)
        self.cell_vol = np.array([0.3] * SIZE)

class Phase(object):
    def __init__(self, q, F, kb, T, tm):
        self.q = q
        self.F = F
        self.T = T
        self.kb = kb
        self.tm = tm
        self.cells = Cells()
        self.p = self
        self.dt = 0.1

p = Phase(0.1, 0.1, 0.1, 0.1, 0.1)

channel_config = [{'apply to': 'all',
  'channel class': 'K',
  'channel inhibitors': ['Anion'],
  'channel type': 'KLeak',
  'inhibitor Km': [0.392],
  'inhibitor n': [3.46],
  'inhibitor zone': ['cell'],
  'max conductivity': 9.15e-17,
  'name': 'K_Channel1'},
 {'apply to': 'all', 'channel class': 'K', 'channel type': 'Kir2p1', 'max Dm': 7.27e-17, 'name': 'Kir2p1'},
 {'apply to': 'all', 'channel class': 'Fun', 'channel type': 'HCN4', 'max Dm': 9.849e-17, 'name': 'HCN4'},
 {'apply to': 'all', 'channel class': 'Fun', 'channel type': 'HCN2', 'max Dm': 7.9e-19, 'name': 'HCN2'},
 {'apply to': 'all', 'channel class': 'Na', 'channel type': 'Nav1p3', 'init active': False, 'max Dm': 1e-14, 'name': 'Nav'},
 {'apply to': 'all', 'channel class': 'K', 'channel type': 'Kv1p5', 'max Dm': 8.16311e-15, 'name': 'Kv'}]

def get_conductivity(D, z, c, d, p):
    """
    Simple function returns conductivity (S/m)
    for diffusion constant D (m2/s), mean
    concentration c, and membrane thickness 'd'.

    """
    # (D * (z ** 2) * p.F * c) / (d * p.R * p.T)
    return (D * p.q *(z**2) * p.F * c) / (d * p.kb * p.T)


class Sim(object):
    def __init__(self, vm):
        self.mdl = SIZE
        self.edl = SIZE
        self.vm = vm
        self.mem_concs = {'Anion': np.ones(SIZE)}
        self.zmol = {
            'K': 0.01,
            'Na': 0.02,
            'Ca': 0.03
        }
        self.iNa = 0.1
        self.iK = 0.1
        self.iCa = 0.1

        self.cbar_dic = {
            'K': np.array([11.25] * SIZE),
            'Na': np.array([62.75] * SIZE),
            'Ca': np.array([13.01] * SIZE)
        }

        self.geo_conv = 0.1

        self.rev_E_dic = {
            'K': np.array([-3.1] * SIZE),
            'Na': np.array([0.07] * SIZE),
            'Ca': np.array([-1.2] * SIZE)
        }
        self.extra_J_mem = np.zeros(SIZE)

    def get_ion(self, ion_name: str) -> int:
        '''
        0-based index assigned to the ion with the passed name (e.g., ``Na``)
        if this ion is enabled by this simulation *or* the empty list
        otherwise.

        This index is guaranteed to uniquely (but arbitrarily) identify this
        ion with respect to this simulation.
        '''

        if ion_name == 'Na':
            ion = self.iNa
        elif ion_name == 'K':
            ion = self.iK
        elif ion_name == 'Ca':
            ion = self.iCa
        elif ion_name == 'Cl':
            ion = self.iCl
        elif ion_name == 'M':
            ion = self.iM
        elif ion_name == 'H':
            ion = self.iH
        elif ion_name == 'P':
            ion = self.iP
        else:
            #FIXME: Shouldn't a fatal exception be raised instead?
            logs.warning(
                'Oops! Molecule gated ion "%s" channel target not found!',
                ion_name)
            ion = []

        return ion

    def run_fast_loop_channels(self, phase):
        '''
        Perform all dynamic channels enabled by this gene regulatory network (GRN)
        for the current simulation time step.

        Parameters
        --------
        phase : SimPhase
            Current simulation phase.
        '''

        sim = phase.sim
        cells = phase.cells
        p = phase.p

        # get the object corresponding to the specific channel:
        for i, name in enumerate(self.channels):

            chan = self.channels[name]

            # compute the channel activity
            # calculate the value of the channel modulation constant:
            moddy = chan.alpha_eval(self, sim, cells, p)

            # set the modulator state in the channel core
            chan.channel_core.modulator = moddy

            # run the channel to update its state:
            chan.channel_core.run(sim.vm, p)

            # update concentrations according to the channel state:
            for ion, rel_perm in zip(chan.channel_core.ions, chan.channel_core.rel_perm):

                # get the charge of the ion:
                zzz = self.zmol[ion]

                ioni = sim.get_ion(ion)

                # obtain the channel's effective diffusion constant:
                DChan = chan.channel_core.P*rel_perm*chan.maxDm*moddy

                G_Chan = get_conductivity(DChan, zzz, sim.cbar_dic[ion], p.tm, p)*sim.geo_conv

                # G_Chan = np.dot(cells.M_sum_mems, G_Chano)/cells.num_mems

                # J_ED = G_Chan * (sim.vm - sim.rev_E_dic[ion]) * (1 / cells.num_mems[cells.mem_to_cells])
                J_ED = G_Chan*(sim.vm - sim.rev_E_dic[ion])

                # Add channel flux to the membrane fluxes data array:
                self.extra_J_mem += J_ED

                # Store channel flux specific to the channel as well:
                chan.channel_core.chan_flux = J_ED/(zzz*p.F)

phase = Phase(0.1, 0.1, 0.1, 0.01, 0.1)
vm = np.ones(SIZE)
sim = Sim(vm)
phase.sim = sim
m = MasterOfNetworks(sim, phase.cells, {}, phase)
m.read_channels(channel_config, phase)
sim.channels = m.channels

SIM_SIZE = 1000

t0 = time.time()
for i in range(SIM_SIZE):
    sim.run_fast_loop_channels(phase)
print("Betse: %.2f" % (time.time() - t0))

def set_relperm(ch, val):
    for i, v in enumerate(val):
        ch.rel_perm[i] = v;

def arr_addr(a):
    return ffi.from_buffer("double*", a)

keepalive = set() # bag of shit that we like to be kept alive and not randomly freed

KLeak = ffi.new("struct Channel*")
KLeak.init = ffi.addressof(lib, 'init_KLeak_channel')
KLeak.run = ffi.addressof(lib, 'run_KLeak_channel')
KLeak.maxDm = 0.1
set_relperm(KLeak, [0, 1.0, 0, 0])
KLeak.is_identity = 0
KLeak.calculate_moddy = ffi.addressof(lib, 'calc_moddy_anion')

Kir2p1 = ffi.new("struct Channel*")
Kir2p1.m = ffi.cast("double*", lib.malloc(8 * SIZE))
Kir2p1.h = ffi.cast("double*", lib.malloc(8 * SIZE))
lib.init_Kir2p1_channel(Kir2p1, arr_addr(vm))
Kir2p1.run = ffi.addressof(lib, 'run_Kir2p1_channel')
Kir2p1.maxDm = 0.1
set_relperm(Kir2p1, [0, 1.0, 0, 0])
Kir2p1.is_identity = 1
Kir2p1.calculate_moddy = ffi.addressof(lib, 'calc_moddy_identity')

HCN4 = ffi.new("struct Channel*")
HCN4.m = ffi.cast("double*", lib.malloc(8 * SIZE))
HCN4.h = ffi.cast("double*", lib.malloc(8 * SIZE))
lib.init_HCN4_channel(HCN4, arr_addr(vm))
HCN4.run = ffi.addressof(lib, 'run_HCN4_channel')
HCN4.maxDm = 0.1
set_relperm(HCN4, [0, 0.2, 1.0, 0.05])
HCN4.is_identity = 1
HCN4.calculate_moddy = ffi.addressof(lib, 'calc_moddy_identity')

HCN2 = ffi.new("struct Channel*")
HCN2.m = ffi.cast("double*", lib.malloc(8 * SIZE))
HCN2.h = ffi.cast("double*", lib.malloc(8 * SIZE))
lib.init_HCN2_channel(HCN2, arr_addr(vm))
HCN2.run = ffi.addressof(lib, 'run_HCN2_channel')
HCN2.maxDm = 0.1
set_relperm(HCN2, [0, 0.2, 1.0, 0.05])
HCN2.is_identity = 1
HCN2.calculate_moddy = ffi.addressof(lib, 'calc_moddy_identity')

Kv1p5 = ffi.new("struct Channel*")
Kv1p5.m = ffi.cast("double*", lib.malloc(8 * SIZE))
Kv1p5.h = ffi.cast("double*", lib.malloc(8 * SIZE))
lib.init_Kv1p5_channel(Kv1p5, arr_addr(vm))
Kv1p5.run = ffi.addressof(lib, 'run_Kv1p5_channel')
Kv1p5.maxDm = 0.1
set_relperm(Kv1p5, [0, 0.2, 0.0, 0.0])
Kv1p5.is_identity = 1
Kv1p5.calculate_moddy = ffi.addressof(lib, 'calc_moddy_identity')

sim = ffi.new("struct Sim*")
p = ffi.new("struct Phase*")
p.q = 0.1
p.F = 0.1
p.kb = 0.1
p.T = 0.01
p.dt = 0.1
p.tm = 0.1
sim.p = p
for i in range(4):
    sim.zmol[i] = 0.1
cbar = [np.zeros(SIZE) for i in range(4)]
for i in range(4):
    sim.cbar_dic[i] = arr_addr(cbar[i])
rev_e = [np.zeros(SIZE) for i in range(4)]
for i in range(4):
    sim.rev_e_dic[i] = arr_addr(rev_e[i])

channels = ffi.new("struct Channel*[]", 5)
channels[0] = KLeak
channels[1] = Kir2p1
channels[2] = HCN4
channels[3] = HCN2
channels[4] = Kv1p5

t0 = time.time()
for i in range(SIM_SIZE):
    lib.run_fast_loop_channels(arr_addr(vm), channels, 5, sim, arr_addr(vm))
print("Fast C: %.2f" % (time.time() - t0))
#KLeak.
if 0:
    p_f = ffi.new("struct Phase*")
    p_f.q = 0.1
    p_f.F = 0.1
    p_f.kb = 0.1
    p_f.T = 0.01


    def get_conductivity(D, z, c, d, p):
        """
        Simple function returns conductivity (S/m)
        for diffusion constant D (m2/s), mean
        concentration c, and membrane thickness 'd'.

        """
        # (D * (z ** 2) * p.F * c) / (d * p.R * p.T)
        return (D * p.q *(z**2) * p.F * c) / (d * p.kb * p.T)

    def arr_addr(a):
        return ffi.from_buffer("double*", a)

    SIZE = 16000

    print("SIZE: %s" % SIZE)
    P = np.zeros(SIZE, dtype=float)
    moddy = np.zeros(SIZE, dtype=float)
    c = np.zeros(SIZE, dtype=float)
    out = np.zeros(SIZE, dtype=float)
    t0 = time.time()
    for i in range(100000):
        run_fast_loop_channels(P, moddy, 0.1, 0.1, c, p)
    dt = time.time() - t0
    print("Numpy: %.2f" % dt)
    t0 = time.time()
    addr_out = arr_addr(out)
    addr_P = arr_addr(P)
    addr_c = arr_addr(c)
    addr_moddy = arr_addr(moddy)
    t0 = time.time()
    for i in range(100000):
        lib.run_fast_loop_channels(addr_out, addr_P, SIZE, addr_moddy, SIZE, 0.1, 0.1, addr_c, SIZE, p_f)
    dt = time.time() - t0
    print("Pure C: %.2f" % (dt,))
