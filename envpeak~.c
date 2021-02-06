/* ---------------- envpeak~ - simple envpeakelope follower. ----------------- */


#include "m_pd.h"
#include "math.h"

#define MAXOVERLAP 32
#define INITVSTAKEN 64

typedef struct sigenvpeak
{
    t_object x_obj;                 /* header */
    void *x_outlet;                 /* a "float" outlet */
    t_outlet *x_outlet_peak;                 /* a "float" outlet */
    void *x_clock;                  /* a "clock" object */
    t_sample *x_buf;                   /* a Hanning window */
    int x_phase;                    /* number of points since last output */
    int x_period;                   /* requested period of output */
    int x_realperiod;               /* period rounded up to vecsize multiple */
    int x_npoints;                  /* analysis window size in samples */
    t_float x_result;                 /* result to output */
    t_float x_result_peak;                 /* peak result */
    t_float x_peakp;                 /* peak of period */
    t_sample x_sumbuf[MAXOVERLAP];     /* summing buffer */
    t_float x_f;
    int x_allocforvs;               /* extra buffer for DSP vector size */
} t_sigenvpeak;

t_class *envpeak_tilde_class;
static void envpeak_tilde_tick(t_sigenvpeak *x);

static void *envpeak_tilde_new(t_floatarg fnpoints, t_floatarg fperiod)
{
    int npoints = fnpoints;
    int period = fperiod;
    t_sigenvpeak *x;
    t_sample *buf;
    int i;

    if (npoints < 1) npoints = 1024;
    if (period < 1) period = npoints/2;
    if (period < npoints / MAXOVERLAP + 1)
        period = npoints / MAXOVERLAP + 1;
    if (!(buf = getbytes(sizeof(t_sample) * (npoints + INITVSTAKEN))))
    {
        error("envpeak: couldn't allocate buffer");
        return (0);
    }
    x = (t_sigenvpeak *)pd_new(envpeak_tilde_class);
    x->x_buf = buf;
    x->x_npoints = npoints;
    x->x_phase = 0;
    x->x_period = period;
    for (i = 0; i < MAXOVERLAP; i++) x->x_sumbuf[i] = 0;
    for (i = 0; i < npoints; i++)
        buf[i] = (1. - cos((2 * 3.14159 * i) / npoints))/npoints;
    for (; i < npoints+INITVSTAKEN; i++) buf[i] = 0;
    x->x_clock = clock_new(x, (t_method)envpeak_tilde_tick);
    x->x_outlet = outlet_new(&x->x_obj, gensym("float"));
    x->x_outlet_peak = outlet_new(&x->x_obj, gensym("float"));
    x->x_f = 0;
    x->x_peakp = 0;
    x->x_allocforvs = INITVSTAKEN;
    return (x);
}

static t_int *envpeak_tilde_perform(t_int *w)
{
    t_sigenvpeak *x = (t_sigenvpeak *)(w[1]);
    t_sample *in = (t_sample *)(w[2]);
    int n = (int)(w[3]);
    int count;
    t_sample *sump;
    in += n;
    for (count = x->x_phase, sump = x->x_sumbuf;
        count < x->x_npoints; count += x->x_realperiod, sump++)
    {
        t_sample *hp = x->x_buf + count;
        t_sample *fp = in;
        t_sample sum = *sump;
        int i;

        for (i = 0; i < n; i++)
        {
            fp--;
            sum += *hp++ * (*fp * *fp);

	    if (*fp > x->x_peakp) x->x_peakp = *fp;
        }
        *sump = sum;
    }
    sump[0] = 0;
    x->x_phase -= n;
    if (x->x_phase < 0)
    {
        x->x_result = x->x_sumbuf[0];
	x->x_result_peak = x->x_peakp;
	x->x_peakp = 0;
        for (count = x->x_realperiod, sump = x->x_sumbuf;
            count < x->x_npoints; count += x->x_realperiod, sump++)
                sump[0] = sump[1];
        sump[0] = 0;
        x->x_phase = x->x_realperiod - n;
        clock_delay(x->x_clock, 0L);
    }
    return (w+4);
}

static void envpeak_tilde_dsp(t_sigenvpeak *x, t_signal **sp)
{
    if (x->x_period % sp[0]->s_n) x->x_realperiod =
        x->x_period + sp[0]->s_n - (x->x_period % sp[0]->s_n);
    else x->x_realperiod = x->x_period;
    if (sp[0]->s_n > x->x_allocforvs)
    {
        void *xx = resizebytes(x->x_buf,
            (x->x_npoints + x->x_allocforvs) * sizeof(t_sample),
            (x->x_npoints + sp[0]->s_n) * sizeof(t_sample));
        if (!xx)
        {
            error("envpeak~: out of memory");
            return;
        }
        x->x_buf = (t_sample *)xx;
        x->x_allocforvs = sp[0]->s_n;
    }
    dsp_add(envpeak_tilde_perform, 3, x, sp[0]->s_vec, (t_int)sp[0]->s_n);
}

static void envpeak_tilde_tick(t_sigenvpeak *x) /* callback function for the clock */
{
    outlet_float(x->x_outlet, powtodb(x->x_result));
    outlet_float(x->x_outlet_peak, x->x_result_peak);
}

static void envpeak_tilde_ff(t_sigenvpeak *x)           /* cleanup on free */
{
    clock_free(x->x_clock);
    freebytes(x->x_buf, (x->x_npoints + x->x_allocforvs) * sizeof(*x->x_buf));
}


void envpeak_tilde_setup(void)
{
    envpeak_tilde_class = class_new(gensym("envpeak~"), (t_newmethod)envpeak_tilde_new,
        (t_method)envpeak_tilde_ff, sizeof(t_sigenvpeak), 0, A_DEFFLOAT, A_DEFFLOAT, 0);
    CLASS_MAINSIGNALIN(envpeak_tilde_class, t_sigenvpeak, x_f);
    class_addmethod(envpeak_tilde_class, (t_method)envpeak_tilde_dsp,
        gensym("dsp"), A_CANT, 0);
}

