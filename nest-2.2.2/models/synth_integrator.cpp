/*
 *  synth_integrator.cpp
 *
 *  This file is part of NEST
 *
 *  Copyright (C) 2004-2008 by
 *  The NEST Initiative
 *
 *  See the file AUTHORS for details.
 *
 *  Permission is granted to compile and modify
 *  this file for non-commercial use.
 *  See the file LICENSE for details.
 *
 *  (C) 2010 Jan Morén 
 */

//#include "scmodules_names.h"
#include "exceptions.h"
#include "synth_integrator.h"
#include "network.h"
#include "dict.h"
#include "integerdatum.h"
#include "doubledatum.h"
#include "dictutils.h"
#include "numerics.h"
//#include "analog_data_logger_impl.h"
#include "universal_data_logger_impl.h"

#include <limits>

/* ---------------------------------------------------------------- 
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap<nest::synth_integrator> nest::synth_integrator::recordablesMap_;

using namespace nest;

namespace nest {
/*
   * Override the create() method with one call to RecordablesMap::insert_() 
   * for each quantity to be recorded.
   */
    template <>
    void RecordablesMap<nest::synth_integrator>::create()
    {
	// use standard names whereever you can for consistency!
	//insert_(nest::names::V_m, &nest::synth_integrator::get_V_m_);
	insert_("weighted_spikes_ex", &nest::synth_integrator::get_weighted_spikes_ex_);
	insert_("weighted_spikes_in", &nest::synth_integrator::get_weighted_spikes_in_);
    }
}
/* ---------------------------------------------------------------- 
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::synth_integrator::Parameters_::Parameters_()
    : Smax_	    ( 1000.0    ),	// integrated value for 1/ms spike
    Gamma_	    (1.0    ),		//  gamma scale
    Reset_	    (10.0    )		//  reset input threshold
{}

nest::synth_integrator::State_::State_()
    : acc_   (0.0)
{}

/* ---------------------------------------------------------------- 
 * Paramater and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void nest::synth_integrator::Parameters_::get(DictionaryDatum &d) const
{

    def<double>(d, names::Smax, Smax_);
    def<double>(d, names::Gamma, Gamma_);
    def<double>(d, names::Reset, Reset_);
}

void nest::synth_integrator::Parameters_::set(const DictionaryDatum& d)
{
    updateValue<double>(d, names::Gamma, Gamma_);
    updateValue<double>(d, names::Reset, Reset_);
    updateValue<double>(d, names::Smax, Smax_);
}

void nest::synth_integrator::State_::get(DictionaryDatum &d, const Parameters_& p) const
{
    //  def<double>(d, nest::names::V_m, y3_ + p.U0_); // Membrane potential
}

void nest::synth_integrator::State_::set(const DictionaryDatum& d, const Parameters_& p)
{
    //if ( updateValue<double>(d, nest::names::V_m, y3_) )
    //  y3_ -= p.U0_;
}
nest::synth_integrator::Buffers_::Buffers_(synth_integrator& n)
  : logger_(n)
{}

nest::synth_integrator::Buffers_::Buffers_(const Buffers_&, synth_integrator& n)
  : logger_(n)
{}
/* ---------------------------------------------------------------- 
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

    nest::synth_integrator::synth_integrator()
: Archiving_Node(), 
    P_(), 
    S_(),
    B_(*this)
{
  recordablesMap_.create();
}

    nest::synth_integrator::synth_integrator(const synth_integrator& n)
: Archiving_Node(n), 
    P_(n.P_), 
    S_(n.S_),
    B_(n.B_, *this)
{}

/* ---------------------------------------------------------------- 
 * Node initialization functions
 * ---------------------------------------------------------------- */

void nest::synth_integrator::init_node_(const Node& proto)
{
    const synth_integrator& pr = downcast<synth_integrator>(proto);
    P_ = pr.P_;
    S_ = pr.S_;
}

void nest::synth_integrator::init_state_(const Node& proto)
{
    const synth_integrator& pr = downcast<synth_integrator>(proto);
    S_ = pr.S_;
}

void nest::synth_integrator::init_buffers_()
{
    B_.ex_spikes_.clear();       // includes resize
    B_.in_spikes_.clear();       // includes resize
    B_.logger_.reset();
    Archiving_Node::clear_history();
}

void nest::synth_integrator::calibrate()
{
    B_.logger_.init();  // ensures initialization in case mm connected after Simulate
    V_.t_res = Time::get_resolution().get_ms(); 
}

/* ---------------------------------------------------------------- 
 * Update and spike handling functions
 */

void nest::synth_integrator::update(Time const & origin, const long_t from, const long_t to)
{
    assert(to >= 0 && (delay) from < Scheduler::get_min_delay());
    assert(from < to);

    double rn, prop;

    double resetc = 0.0;

    for ( long_t lag = from ; lag < to ; ++lag )
    {

	V_.weighted_spikes_ex_ = B_.ex_spikes_.get_value(lag);	
	V_.weighted_spikes_in_ = B_.in_spikes_.get_value(lag);

	S_.acc_ += V_.weighted_spikes_ex_;
	

	resetc += (-V_.weighted_spikes_in_);
	if (resetc > P_.Reset_)
	    S_.acc_ = 0.0;

	prop = pow((S_.acc_/(P_.Smax_)),(1.0/P_.Gamma_));
	rn = random()/(RAND_MAX+1.0)/V_.t_res;

	if (rn<= prop) { 
	    SpikeEvent se;
	    network()->send(*this, se, lag);
	}
    
	B_.logger_.record_data(origin.get_steps() + lag);
    }                           
}


void nest::synth_integrator::handle(SpikeEvent & e)
{
    assert(e.get_delay() > 0);

    //printf ("\tSynapse: %d\n", e.get_rport());
    if (e.get_rport()==0) {
	if(e.get_weight() > 0.0)
	    B_.ex_spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
		    e.get_weight() * e.get_multiplicity() );
	else if (e.get_weight() < 0.0)
	    B_.in_spikes_.add_value(e.get_rel_delivery_steps(network()->get_slice_origin()),
		    e.get_weight() * e.get_multiplicity() );
    }
}

void nest::synth_integrator::handle(CurrentEvent& e)
{
    assert(e.get_delay() > 0);
}

void nest::synth_integrator::handle(DataLoggingRequest& e)
{
    B_.logger_.handle(e);
}

// namespace
