typedef struct Phase {
	double kb;
	double q;
	double F;
	double T;
	double tm;
	double dt;
} Phase;

#define SIZE 16482

#define ION_NONE 0
#define ION_K 1
#define ION_NA 2
#define ION_CA 3
#define ION_LAST 4

struct Channel;

typedef struct Channel {
	int(*run)(struct Channel*, double *, double *, struct Phase*);
	int(*init)(struct Channel*, double *);
	double(*calculate_moddy)(double, int);
	int is_identity;
	double maxDm;
	double rel_perm[ION_LAST];

	double *m;
	double *h;
} Channel;

typedef struct Sim {
	double zmol[ION_LAST];
	double* cbar_dic[ION_LAST];
	double* rev_e_dic[ION_LAST];
	double geo_conv;
	struct Phase* p;
} Sim;

extern double Anion[SIZE];

int get_conductivity(double *out, double *D, double z, double* c, double d, struct Phase*, double geo_conv);

int run_fast_loop_channels(double *out, struct Channel* channels[], int num_channels, struct Sim* sim, double *vm);

// specific channel functions
int init_KLeak_channel(struct Channel *chan, double *vm);
int run_KLeak_channel(struct Channel *chan, double *P, double *vm, struct Phase *p);
int init_Kir2p1_channel(struct Channel *chan, double *vm);
int run_Kir2p1_channel(struct Channel* chan, double *P, double *vm, struct Phase *p);
int init_HCN4_channel(struct Channel *chan, double *vm);
int run_HCN4_channel(struct Channel* chan, double *P, double *vm, struct Phase *p);
int init_HCN2_channel(struct Channel* chan, double *vm);
int run_HCN2_channel(struct Channel *chan, double *P, double *vm, struct Phase *p);
int init_Kv1p5_channel(struct Channel* chan, double *vm);
int run_Kv1p5_channel(struct Channel *chan, double *P, double *vm, struct Phase *p);

// specific moddy functions
double calc_moddy_identity(double v, int i);
double calc_moddy_anion(double v, int i);

