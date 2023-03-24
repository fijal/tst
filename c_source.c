
#include "sim_tb.h"

#include <math.h>
#include <assert.h>
#include <stdlib.h>
#include <stdio.h>

double Anion[SIZE];

int init_KLeak_channel(struct Channel *chan, double *vm)
{
	return 0;
}

int run_KLeak_channel(struct Channel *chan, double *P, double *vm, struct Phase *p)
{
    for (int i = 0; i < SIZE; ++i) {
    	P[i] = 1.0;
    }
	return 0;
}

int init_Kir2p1_channel(struct Channel *chan, double *vm)
{
    // initialize values of the m and h gates of the sodium channel based on m_inf and h_inf:
    for (int i = 0; i < SIZE; i++) {
	    chan->m[i] = 1 / (1 + exp((vm[i] - (-96.48)) / 23.26));
	    chan->h[i] = 1 / (1 + exp((vm[i] - (-168.28)) / -44.13));
	}
	return 0;
}

int run_Kir2p1_channel(struct Channel* chan, double *P, double *vm, struct Phase *p)
{
	double time_unit = 1.0e3;

	for (int i = 0; i < SIZE; i++) {

		double _mInf = 1 / (1 + exp((vm[i] * 1000 - (-96.48)) / 23.26));
        double _mTau = 3.7 + (-3.37 / (1 + exp((vm[i] - -32.9) / 27.93)));
        double _hInf = 1 / (1 + exp((vm[i] * 1000 - (-168.28)) / -44.13));
       	double _hTau = 0.85 + (306.3 / (1 + exp((vm[i] - -118.29) / -27.23)));

        double dt = p->dt*time_unit;

        chan->m[i] = (_mTau * chan->m[i] + dt * _mInf)/(_mTau + dt);
        chan->h[i] = (_hTau * chan->h[i] + dt * _hInf)/(_hTau + dt);


        P[i] = chan->m[i] + (chan->h[i] * chan->h[i]);
	}
	return 0;
}

int init_HCN4_channel(struct Channel *chan, double *vm)
{
    for (int i = 0; i < SIZE; i++) {
    	chan->m[i] = 1.000 / (1 + exp(vm[i] - 10 - (-100))/9.6);
    }
    return 0;
}

int run_HCN4_channel(struct Channel* chan, double *P, double *vm, struct Phase *p)
{
    double dt = p->dt * 1e03;
    double _mTau = 461.0;

    for (int i = 0; i < SIZE; i++) {
    	double _mInf = 1.000 / (1 + exp(vm[i] * 1000 - -100) / 9.6);
        P[i] = (_mTau*chan->m[i] + dt*_mInf)/(_mTau + dt);
    }
    return 0;
}

int init_HCN2_channel(struct Channel* chan, double *vm)
{
    for (int i = 0; i < SIZE; i++) {
    	chan->m[i] = 1.000 / (1 + exp(vm[i] - 10 - (-99))/6.2);
    }
	return 0;
}

int run_HCN2_channel(struct Channel *chan, double *P, double *vm, struct Phase *p)
{
    double dt = p->dt * 1e03;
    double _mTau = 184.0;

    for (int i = 0; i < SIZE; i++) {
    	double _mInf = 1.000 / (1 + exp(vm[i] * 1000 - -99) / 6.2);
        P[i] = (_mTau*chan->m[i] + dt*_mInf)/(_mTau + dt);
    }
    return 0;
}

int init_Kv1p5_channel(struct Channel* chan, double *vm)
{
    for (int i = 0; i < SIZE; i++) {
        chan->m[i] = 1.0000 / (1 + exp((vm[i] - -6.0000) / -6.4000));
        chan->h[i] = 1.0000 / (1 + exp((vm[i] - -25.3000) / 3.5000));
    }
    return 0;
}

int run_Kv1p5_channel(struct Channel *chan, double *P, double *vm, struct Phase *p)
{
	double time_unit = 1e3;

	for (int i = 0; i < SIZE; i++) {
        double _mInf = 1.0000 / (1 + exp((vm[i] - -6.0000) / -6.4000));
        double _mTau = (-0.1163 * vm[i]) + 8.3300;
        double _hInf = 1.0000 / (1 + exp((vm[i] - -25.3000) / 3.5000));
        double _hTau = (-15.5000 * vm[i]) + 1620.0000;

        double dt = p->dt*time_unit;

        chan->m[i] = (_mTau * chan->m[i] + dt * _mInf) / (_mTau + dt);
        chan->h[i] = (_hTau * chan->h[i] + dt * _hInf) / (_hTau + dt);

        P[i] = chan->m[i] + chan->h[i];
    }
    return 0;
}

double calc_moddy_identity(double v, int i)
{
	return v;
}

double calc_moddy_anion(double v, int i)
{
	return v * (1 / (1 + pow(Anion[i] / 0.392, 3.46)));
}

int get_conductivity(double *out, double *D, double z, double* c, double d, struct Phase* p, double geo_conv)
{
    // return (D * p.q *(z**2) * p.F * c) / (d * p.kb * p.T)
	for (int i = 0; i < SIZE; i++) {
		out[i] = ((D[i] * p->q * (z * z) * p->F * c[i]) / d * p->kb * p->T) * geo_conv;
	}
	return 0;
}

int compute_dchan(double *DChan, double *P, double rel_perm, double maxDm, int is_identity, double(*calc_moddy)(double, int))
{
	for (int i = 0; i < SIZE; ++i)
	{
		if (is_identity) {
			DChan[i] = P[i] * rel_perm * maxDm;
		} else {
			DChan[i] = calc_moddy_anion(P[i] * rel_perm * maxDm, i);
		}
	}
	return 0;
}

int compute_JED(double *out, double *GChan, double*vm, double* rev_e)
{
	for (int i = 0; i < SIZE; ++i)
	{
		out[i] = GChan[i] * (vm[i] - rev_e[i]);
	}
	return 0;
}

int run_fast_loop_channels(double *out, struct Channel* channels[], int num_channels, struct Sim* sim, double *vm)
{
	double P[SIZE];
	double DChan[SIZE];

	for (int chan_no = 0; chan_no < num_channels; chan_no++) {
		struct Channel* chan = channels[chan_no];

		chan->run(chan, P, vm, sim->p);
		for (int ion_no = 0; ion_no < ION_LAST; ++ion_no) {
			if (chan->rel_perm[ion_no] == 0)
				continue;
			double zzz = sim->zmol[ion_no];
			compute_dchan(P, P, chan->rel_perm[ion_no], chan->maxDm, chan->is_identity, chan->calculate_moddy);
			get_conductivity(DChan, DChan, zzz, sim->cbar_dic[ion_no], sim->p->tm, sim->p, sim->geo_conv);
			compute_JED(out, DChan, vm, sim->rev_e_dic[ion_no]);
		}
	}
	return 0;
}