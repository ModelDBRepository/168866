/*
 *  aeif_cond_nmda_alpha.cpp
 *
 * Adaptive exponential integrate-and-fire model (Brette and Gerstner) with
 * alpha-shaped synaptic conductances, as implemented in aeif_cond_alpha.cpp 
 * from the nest simulator (v. 1.9.8498).
 *
 * This version adds an NMDA synapse-like input with positive
 * voltage-dependence on the conductance.
 *
 * Jan Moren, 2009
 *
 * 2011 - Update for nest2
 *
 * Original model:
*
 *  This file is part of NEST
 *
 *  Copyright (C) 2005 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 */


#include "aeif_cond_nmda_alpha.h"
//#include "scmodules_names.h"
#include "nest_names.h"

#ifdef HAVE_GSL_1_11

#include "universal_data_logger_impl.h"
#include "exceptions.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
#include <limits>

#include <cmath>
#include <iomanip>
#include <iostream>
#include <cstdio>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<nest::aeif_cond_nmda_alpha> nest::aeif_cond_nmda_alpha::recordablesMap_;

using namespace nest;

namespace nest {
/*
   * Override the create() method with one call to RecordablesMap::insert_() 
   * for each quantity to be recorded.
   */
    template <>
    void RecordablesMap<aeif_cond_nmda_alpha>::create()
    {
	// use standard names whereever you can for consistency!
	insert_(names::V_m, &aeif_cond_nmda_alpha::get_y_elem_<aeif_cond_nmda_alpha::V_M>);
	insert_(names::g_ex, &aeif_cond_nmda_alpha::get_y_elem_<aeif_cond_nmda_alpha::G_EX>);
	insert_(names::g_in, &aeif_cond_nmda_alpha::get_y_elem_<aeif_cond_nmda_alpha::G_IN>);
	insert_(names::g_n, &aeif_cond_nmda_alpha::get_y_elem_<aeif_cond_nmda_alpha::G_N>);
	insert_(names::w, &aeif_cond_nmda_alpha::get_y_elem_<aeif_cond_nmda_alpha::W>);
	//insert_("weighted_spikes_ex", &nest::synth_integrator::get_weighted_spikes_ex_);
	//insert_("weighted_spikes_in", &nest::synth_integrator::get_weighted_spikes_in_);
    }
}

 extern "C"
 int aeif_cond_nmda_alpha_dynamics (double, const double y[], double f[], void* param)
 {
   // shorthand for class we work for
   typedef nest::aeif_cond_nmda_alpha AEIF;
   
   // get parameters as reference  
   AEIF::Parameters_* tmp =
     reinterpret_cast<AEIF::Parameters_*>(param);
   assert(tmp);
   AEIF::Parameters_& p = *tmp;

   // shorthand for state variables
   const nest::double_t& V     = y[AEIF::V_M  ];
   const nest::double_t& dg_ex = y[AEIF::DG_EX];
   const nest::double_t&  g_ex = y[AEIF::G_EX ];
   const nest::double_t& dg_in = y[AEIF::DG_IN];
   const nest::double_t&  g_in = y[AEIF::G_IN ];
   const nest::double_t& w     = y[AEIF::W    ];
   const nest::double_t& dg_n = y[AEIF::DG_N];
   const nest::double_t&  g_n = y[AEIF::G_N ];

   const nest::double_t I_syn_exc = g_ex * (V - p.E_ex);
   const nest::double_t I_syn_inh = g_in * (V - p.E_in);
   const nest::double_t I_syn_n   = g_n * (V - p.E_n);
   const nest::double_t I_spike = p.Delta_T * std::exp((V - p.V_th) / p.Delta_T);

   // dv/dt
   f[AEIF::V_M  ] = ( -p.g_L *( (V-p.E_L) - I_spike) 
 	                   - I_syn_exc - I_syn_inh - I_syn_n - w + p.I_e + p.I_stim) / p.C_m;

   f[AEIF::DG_EX] = -dg_ex / p.tau_syn_ex;
   f[AEIF::G_EX ] =  dg_ex - g_ex / p.tau_syn_ex; // Synaptic Conductance (nS)

   f[AEIF::DG_IN] = -dg_in / p.tau_syn_in;
   f[AEIF::G_IN ] =  dg_in - g_in / p.tau_syn_in; // Synaptic Conductance (nS)

   f[AEIF::DG_N] = -dg_n / p.tau_syn_n;
   f[AEIF::G_N ] =  dg_n - g_n / p.tau_syn_n; // Synaptic Conductance (nS)

   // Adaptation current w.
   f[AEIF::W    ] = ( p.a * (V - p.E_L) - w ) / p.tau_w;

   return GSL_SUCCESS;
 }

/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */
    
nest::aeif_cond_nmda_alpha::Parameters_::Parameters_()
  : V_peak_    (  0.0    ),  // mV
    V_reset_   (-60.0    ),  // mV
    t_ref_     (  0.0    ),  // ms
    g_L        ( 30.0    ),  // nS
    C_m        (281.0    ),  // pF
    E_ex       (  0.0    ),  // mV
    E_in       (-85.0    ),  // mV
    E_L        (-70.6    ),  // mV
    Delta_T    (  2.0    ),  // mV
    tau_w      (144.0    ),  // ms
    a          (  4.0    ),  // nS
    b          ( 80.5    ),  // pA
    V_th       (-50.4    ),  // mV
    tau_syn_ex (  0.2    ),  // ms
    tau_syn_in (  2.0    ),  // ms
    I_e        (  0.0    ),  // pA
    
    N_V_min     (-70.0    ),  // mv
    N_V_max     (-50.0    ),  // mV
    N_gain      (  3.0    ),   // -
    tau_syn_n   (  3.0    ),  // ms
    E_n         (  0.0    )  // mV
{}

nest::aeif_cond_nmda_alpha::State_::State_(const Parameters_& p)
  : r_(0)
{
  y_[0] = p.E_L;
  for ( size_t i = 1 ; i < aeif_cond_nmda_alpha::NSTATES ; ++i )
    y_[i] = 0;
}

nest::aeif_cond_nmda_alpha::State_::State_(const State_& s)
  : r_(s.r_)
{
  for ( size_t i = 0 ; i < aeif_cond_nmda_alpha::NSTATES ; ++i )
    y_[i] = s.y_[i];
}

nest::aeif_cond_nmda_alpha::State_& nest::aeif_cond_nmda_alpha::State_::operator=(const State_& s)
{
  assert(this != &s);  // would be bad logical error in program
  
  for ( size_t i = 0 ; i < aeif_cond_nmda_alpha::NSTATES ; ++i )
    y_[i] = s.y_[i];
  r_ = s.r_;
  return *this;
}

/* ---------------------------------------------------------------- 
 * Paramater and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::aeif_cond_nmda_alpha::Parameters_::get(DictionaryDatum &d) const
{
  def<double>(d,names::C_m,    C_m);
  def<double>(d,names::V_th,   V_th);
  def<double>(d,names::t_ref,  t_ref_);
  def<double>(d,names::g_L,    g_L);
  def<double>(d,names::E_L,    E_L); 
  def<double>(d,names::V_reset,V_reset_);
  def<double>(d,names::E_ex,   E_ex);
  def<double>(d,names::E_in,   E_in);
  def<double>(d,names::tau_syn_ex, tau_syn_ex);
  def<double>(d,names::tau_syn_in, tau_syn_in);
  def<double>(d,names::a,      a);
  def<double>(d,names::b,      b);
  def<double>(d,names::Delta_T,Delta_T);
  def<double>(d,names::tau_w,  tau_w);
  def<double>(d,names::I_e,    I_e);
  def<double>(d,names::V_peak, V_peak_);

  def<double>(d, names::N_V_min, N_V_min);
  def<double>(d, names::N_V_max, N_V_max);
  def<double>(d, names::N_gain, N_gain);
  def<double>(d, names::E_n,   E_n);
  def<double>(d, names::tau_syn_n, tau_syn_n);
}

void nest::aeif_cond_nmda_alpha::Parameters_::set(const DictionaryDatum& d)
{
  updateValue<double>(d,names::V_th,    V_th);
  updateValue<double>(d,names::V_peak,  V_peak_);
  updateValue<double>(d,names::t_ref,   t_ref_);
  updateValue<double>(d,names::E_L,     E_L);
  updateValue<double>(d,names::V_reset, V_reset_);
  updateValue<double>(d,names::E_ex,    E_ex);
  updateValue<double>(d,names::E_in,    E_in);
    
  updateValue<double>(d,names::C_m,     C_m);
  updateValue<double>(d,names::g_L,     g_L);
    
  updateValue<double>(d,names::tau_syn_ex, tau_syn_ex);
  updateValue<double>(d,names::tau_syn_in, tau_syn_in);
    
  updateValue<double>(d,names::a,      a);
  updateValue<double>(d,names::b,      b);
  updateValue<double>(d,names::Delta_T,Delta_T);
  updateValue<double>(d,names::tau_w,  tau_w);

  updateValue<double>(d,names::I_e,    I_e);

  updateValue<double>(d, names::N_V_min, N_V_min);
  updateValue<double>(d, names::N_V_max, N_V_max);
  updateValue<double>(d, names::N_gain, N_gain);
  updateValue<double>(d, names::E_n,    E_n);
  updateValue<double>(d, names::tau_syn_n, tau_syn_n);

  if ( V_reset_ >= V_th )
    throw BadProperty("Reset potential must be smaller than threshold.");
    
  if ( V_peak_ <= V_th )
    throw BadProperty("V_peak must be larger than threshold.");

  if ( C_m <= 0 )
  {
    throw BadProperty("Capacitance must be strictly positive.");
  }
    
  if ( t_ref_ < 0 )
    throw BadProperty("Refractory time cannot be negative.");
      
  if ( tau_syn_ex <= 0 || tau_syn_in <= 0 || tau_syn_n <=0 || tau_w <= 0 )
    throw BadProperty("All time constants must be strictly positive.");
}

void nest::aeif_cond_nmda_alpha::State_::get(DictionaryDatum &d) const
{
  def<double>(d,names::V_m,    y_[V_M]);
  def<double>(d,names::g_ex,   y_[G_EX]);
  def<double>(d,names::dg_ex,  y_[DG_EX]);
  def<double>(d,names::g_in,   y_[G_IN]);
  def<double>(d,names::dg_in,  y_[DG_IN]);
  def<double>(d,names::w,      y_[W]);
  def<double>(d,names::g_n,   y_[G_N]);
  def<double>(d,names::dg_n,  y_[DG_N]);
}

void nest::aeif_cond_nmda_alpha::State_::set(const DictionaryDatum& d, const Parameters_&)
{
  updateValue<double>(d,names::V_m,   y_[V_M]);
  updateValue<double>(d,names::g_ex,  y_[G_EX]);
  updateValue<double>(d,names::dg_ex, y_[DG_EX]);
  updateValue<double>(d,names::g_in,  y_[G_IN]);
  updateValue<double>(d,names::dg_in, y_[DG_IN]);
  updateValue<double>(d,names::w,     y_[W]);

  updateValue<double>(d,names::g_n,  y_[G_N]);
  updateValue<double>(d,names::dg_n, y_[DG_N]);

  if ( y_[G_EX] < 0 || y_[G_IN] < 0 || y_[G_N] < 0 )
    throw BadProperty("Conductances must not be negative.");
}

nest::aeif_cond_nmda_alpha::Buffers_::Buffers_(aeif_cond_nmda_alpha& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // The other member variables are left uninitialised or are
  // automatically initialised by their default constructor.
}

nest::aeif_cond_nmda_alpha::Buffers_::Buffers_(const Buffers_&, aeif_cond_nmda_alpha& n)
  : logger_(n),
    s_(0),
    c_(0),
    e_(0)
{
  // The other member variables are left uninitialised or are
  // automatically initialised by their default constructor.
}


/* ---------------------------------------------------------------- 
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::aeif_cond_nmda_alpha::aeif_cond_nmda_alpha()
  : Archiving_Node(), 
    P_(), 
    S_(P_),
    B_(*this)
{
    recordablesMap_.create();
}

nest::aeif_cond_nmda_alpha::aeif_cond_nmda_alpha(const aeif_cond_nmda_alpha& n)
  : Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{
}

nest::aeif_cond_nmda_alpha::~aeif_cond_nmda_alpha()
{
  // GSL structs only allocated by init_nodes_(), so we need to protect destruction
  if ( B_.s_ ) gsl_odeiv_step_free(B_.s_);
  if ( B_.c_ ) gsl_odeiv_control_free(B_.c_);
  if ( B_.e_ ) gsl_odeiv_evolve_free(B_.e_);
}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::aeif_cond_nmda_alpha::init_node_(const Node& proto)
{
  const aeif_cond_nmda_alpha& pr = downcast<aeif_cond_nmda_alpha>(proto);
  P_ = pr.P_;
  S_ = pr.S_;
}

void nest::aeif_cond_nmda_alpha::init_state_(const Node& proto)
{
  const aeif_cond_nmda_alpha& pr = downcast<aeif_cond_nmda_alpha>(proto);
  S_ = pr.S_;
}

void nest::aeif_cond_nmda_alpha::init_buffers_()
{
  B_.spike_exc_.clear();       // includes resize
  B_.spike_inh_.clear();       // includes resize
  B_.spike_nmda_.clear();       // includes resize
  B_.currents_.clear();        // includes resize
//  B_.potentials_.clear_data(); // includes resize
//  B_.conductances_.clear_data(); // includes resize
//  B_.aeif_ws_.clear_data();      // includes resize
  Archiving_Node::clear_history();
  
  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();

  // We must integrate this model with high-precision to obtain decent results
  B_.IntegrationStep_ = std::min(0.01, B_.step_);

  static const gsl_odeiv_step_type* T1 = gsl_odeiv_step_rkf45;
  
  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc (T1, aeif_cond_nmda_alpha::NSTATES);
  else 
    gsl_odeiv_step_reset(B_.s_);
    
  if ( B_.c_ == 0 )  
    B_.c_ = gsl_odeiv_control_yp_new (1e-6,1e-6);
  else
    gsl_odeiv_control_init(B_.c_, 1e-6, 1e-6, 0.0, 1.0);
    
  if ( B_.e_ == 0 )  
    B_.e_ = gsl_odeiv_evolve_alloc(aeif_cond_nmda_alpha::NSTATES);
  else 
    gsl_odeiv_evolve_reset(B_.e_);
  
  B_.sys_.function  = aeif_cond_nmda_alpha_dynamics; 
  B_.sys_.jacobian  = NULL;
  B_.sys_.dimension =  aeif_cond_nmda_alpha::NSTATES;
  B_.sys_.params    = reinterpret_cast<void*>(&P_);
}

void nest::aeif_cond_nmda_alpha::calibrate()
{
      B_.logger_.init();  // ensures initialization in case mm connected after Simulate
  V_.g0_ex_  = 1.0 * numerics::e / P_.tau_syn_ex;
  V_.g0_in_  = 1.0 * numerics::e / P_.tau_syn_in;
  V_.g0_n_  = 1.0 * numerics::e / P_.tau_syn_n;
  V_.RefractoryCounts_ = Time(Time::ms(P_.t_ref_)).get_steps();
  assert(V_.RefractoryCounts_ >= 0);  // since t_ref_ >= 0, this can only fail in error
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void nest::aeif_cond_nmda_alpha::update(Time const & origin, const long_t from, const long_t to)
{
  assert ( to >= 0 && (delay) from < Scheduler::get_min_delay() );
  assert ( from < to );
  assert ( V_M == 0 );

  for ( long_t lag = from; lag < to; ++lag )
  {
    double t = 0.0;

    if ( S_.r_ > 0 )
      --S_.r_;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t; this is of advantage
    // for a consistent and efficient integration across subsequent
    // simulation intervals
    while ( t < B_.step_ )
    {
      const int status = gsl_odeiv_evolve_apply(B_.e_, B_.c_, B_.s_, 
		  	   &B_.sys_,             // system of ODE
			   &t,                   // from t
			    B_.step_,            // to t <= step
			   &B_.IntegrationStep_, // integration step size
			    S_.y_);              // neuronal state

      if ( status != GSL_SUCCESS )
        throw GSLSolverFailure(get_name(), status);

      // spikes are handled inside the while-loop
      // due to spike-driven adaptation
      if ( S_.r_ > 0 )
        S_.y_[V_M] = P_.V_reset_;
      else if ( S_.y_[V_M] >= P_.V_peak_ )
      {
	S_.y_[V_M]  = P_.V_reset_;
	S_.y_[W]   += P_.b; // spike-driven adaptation
	S_.r_       = V_.RefractoryCounts_;
	      
	set_spiketime(Time::step(origin.get_steps() + lag + 1));
	SpikeEvent se;
	network()->send(*this, se, lag);
      }
    }

    S_.y_[DG_EX] += B_.spike_exc_.get_value(lag) * V_.g0_ex_;
    S_.y_[DG_IN] += B_.spike_inh_.get_value(lag) * V_.g0_in_;
     
    // NMDA input

    
    S_.y_[DG_N] += B_.spike_nmda_.get_value(lag) /
	(1 + std::exp(-4.0 * P_.N_gain * 
		      ((S_.y_[V_M] - P_.N_V_min)/(P_.N_V_max - P_.N_V_min) - 0.5))) * 
	V_.g0_n_;


    // set new input current
    P_.I_stim = B_.currents_.get_value(lag);

    // voltage logging
//    B_.potentials_.record_data(origin.get_steps()+lag, S_.y_[V_M]);
//    B_.conductances_.record_data(origin.get_steps()+lag, 
//				   std::pair<nest::double_t, nest::double_t>(S_.y_[G_EX], S_.y_[G_IN]));
//    B_.aeif_ws_.record_data(origin.get_steps()+lag, S_.y_[W]);
    
        B_.logger_.record_data(origin.get_steps() + lag);
  }
}
  
void nest::aeif_cond_nmda_alpha::handle(SpikeEvent & e)
{
  assert(e.get_delay() > 0);
//printf ("\tSynapse: %d\n", e.get_rport());
  if (e.get_rport() == 0) {

      if(e.get_weight() > 0.0)
	B_.spike_exc_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			     e.get_weight() * e.get_multiplicity());
      else
	B_.spike_inh_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			     -e.get_weight() * e.get_multiplicity());  // keep conductances positive

  } else if (e.get_rport() == 1) {
	B_.spike_nmda_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
			     e.get_weight() * e.get_multiplicity());
  }
}

void nest::aeif_cond_nmda_alpha::handle(CurrentEvent& e)
{
  assert(e.get_delay() > 0);

  const nest::double_t c=e.get_current();
  const nest::double_t w=e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()), 
		      w *c);
}

void nest::aeif_cond_nmda_alpha::handle(DataLoggingRequest& e)
{
      B_.logger_.handle(e);
}


#endif //HAVE_GSL_1_11
