/*
 *  synth_integrator.h
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
 *
 * This version is for nest2 prereleases. The logging API has changed so we
 * can no longer use the same module code for different nest versions.
 *
 */

#ifndef SYNTH_INTEGRATOR
#define SYNTH_INTEGRATOR

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"
#include "universal_data_logger.h"
#include "recordables_map.h"

//#include "analog_data_logger.h"

namespace nest
{

  class Network;

  /* BeginDocumentation
Name: synth_integrator - Input integrator.

Description:

  synth integrator takes weighted spiking input, sums it over time, and emits 
  spikes at a rate proportional to the summed input. 

  This is intended as a synthetic replacement or placeholder to a real spike
  integrator circuit. Useful for testing networks where an integrating input 
  is assumed, or to tune a real integrating circuit.

  Hacked-up from iaf_psc_nmda_alpha

Parameters: 

  The following parameters can be set in the status dictionary.

    Smax    double	The accumulated spike*weight at which point the 
			model emits 1 spike/ms.
    Gamma   double	nonlinear scaling factor
    Reset   double	Negative rate input for resetting the integrator


Receives: SpikeEvent, DataLoggingRequest

Author:  February 2010, Moren
*/

  /**
   * Leaky integrate-and-fire neuron with alpha-shaped PSCs.
   */
  class synth_integrator : public Archiving_Node
    {

	public:        

	    typedef Node base;

	    synth_integrator();
	    synth_integrator(const synth_integrator&);

	    /**
	     * Import sets of overloaded virtual functions.
	     * We need to explicitly include sets of overloaded
	     * virtual functions into the current scope.
	     * According to the SUN C++ FAQ, this is the correct
	     * way of doing things, although all other compilers
	     * happily live without.
	     */
/*#ifndef IS_BLUEGENE
	    using Node::check_connection;
#endif*/
	    using Node::connect_sender;
	    using Node::handle;

	    port check_connection(Connection&, port);

	    void handle(SpikeEvent &);
	    void handle(CurrentEvent &);
	    //void handle(PotentialRequest &);
	    void handle(DataLoggingRequest &);

	    port connect_sender(SpikeEvent&, port);
	    port connect_sender(CurrentEvent&, port);
	    //port connect_sender(PotentialRequest&, port);
	    port connect_sender(DataLoggingRequest&, port);

	    void get_status(DictionaryDatum &) const;
	    void set_status(const DictionaryDatum &);

	private:

	    void init_node_(const Node& proto);
	    void init_state_(const Node& proto);
	    void init_buffers_();
	    void calibrate();

	    void update(Time const &, const long_t, const long_t);
	    friend class RecordablesMap<synth_integrator>;
	    friend class UniversalDataLogger<synth_integrator>;

	    // ---------------------------------------------------------------- 

	    /** 
	     * Independent parameters of the model. 
	     */
	    struct Parameters_ {

		/** event limit - at this point we do one spike per ms **/
		double_t Smax_;

		/** gamma factor for probability scaling **/
		double_t Gamma_;
		double_t Reset_;



		Parameters_();  //!< Sets default parameter values

		void get(DictionaryDatum&) const;  //!< Store current values in dictionary
		void set(const DictionaryDatum&);  //!< Set values from dicitonary
	    };

	    // ---------------------------------------------------------------- 

	    /**
	     * State variables of the model.
	     */
	    struct State_ {
		double_t acc_; // Currently accumulated events 

		State_();  //!< Default initialization

		void get(DictionaryDatum&, const Parameters_&) const;
		void set(const DictionaryDatum&, const Parameters_&);
	    };    

	    // ---------------------------------------------------------------- 

	    /**
	     * Buffers of the model.
	     */
	    struct Buffers_ {
		Buffers_(synth_integrator&);
		Buffers_(const Buffers_&, synth_integrator&);
		/** buffers and summs up incoming spikes/currents */
		RingBuffer ex_spikes_;
		RingBuffer in_spikes_;
		//RingBuffer currents_;

		/** Buffer for membrane potential. */
		//AnalogDataLogger<PotentialRequest> potentials_;
		UniversalDataLogger<synth_integrator> logger_;
	    };

	    // ---------------------------------------------------------------- 

	    /**
	     * Internal variables of the model.
	     */
	    struct Variables_ { 

		double_t t_res;	    // Time resolution
		double_t weighted_spikes_ex_; // Not sure we need or want these
		double_t weighted_spikes_in_;

	    };

	    // Access functions for UniversalDataLogger -------------------------------

	    //! Read out the real membrane potential
//	    double_t get_V_m_() const { return S_.y3_ + P_.U0_; }

	    double_t get_weighted_spikes_ex_() const { return V_.weighted_spikes_ex_; }
	    double_t get_weighted_spikes_in_() const { return V_.weighted_spikes_in_; }

	    // ---------------------------------------------------------------- 

	    /**
	     * @defgroup synth_integrator_data
	     * Instances of private data structures for the different types
	     * of data pertaining to the model.
	     * @note The order of definitions is important for speed.
	     * @{
	     */   
	    Parameters_ P_;
	    State_      S_;
	    Variables_  V_;
	    Buffers_    B_;
	    /** @} */
	    
	    static RecordablesMap<synth_integrator> recordablesMap_;
    };

  inline
      port synth_integrator::check_connection(Connection& c, port receptor_type)
      {
	  // printf ("\tCheckSynapse: %d\n", receptor_type);

	  SpikeEvent e;
	  e.set_sender(*this);
	  c.check_event(e);
	  return c.get_target()->connect_sender(e, receptor_type);
      }

  inline
      port synth_integrator::connect_sender(SpikeEvent&, port receptor_type)
      {
	  /* 0: normal receptor; 1: NMDA receptor */
	  //  printf ("\n\tSynapse: %d \n\n", receptor_type);

	  //if (receptor_type >= 1)
	  if (receptor_type != 0)
	      throw UnknownReceptorType(receptor_type, get_name());
	  return receptor_type;
      }

  inline
      port synth_integrator::connect_sender(CurrentEvent&, port receptor_type)
      {
	  if (receptor_type != 0)
	      throw UnknownReceptorType(receptor_type, get_name());
	  return 0;
      }

  inline
      port synth_integrator::connect_sender(DataLoggingRequest& dlr, port receptor_type)
      {
	  if (receptor_type != 0)
	      throw UnknownReceptorType(receptor_type, get_name());
	  //B_.potentials_.connect_logging_device(pr);
	  //return 0;
	  return B_.logger_.connect_logging_device(dlr, recordablesMap_);
      }
  inline
      void synth_integrator::get_status(DictionaryDatum &d) const
      {
	  P_.get(d);
	  S_.get(d, P_);
	  Archiving_Node::get_status(d);
	  (*d)[names::recordables] = recordablesMap_.get_list();
      }

  inline
      void synth_integrator::set_status(const DictionaryDatum &d)
      {
	  Parameters_ ptmp = P_;  // temporary copy in case of errors
	  ptmp.set(d);                       // throws if BadProperty
	  State_      stmp = S_;  // temporary copy in case of errors
	  stmp.set(d, ptmp);                 // throws if BadProperty

	  // We now know that (ptmp, stmp) are consistent. We do not 
	  // write them back to (P_, S_) before we are also sure that 
	  // the properties to be set in the parent class are internally 
	  // consistent.
	  Archiving_Node::set_status(d);

	  // if we get here, temporaries contain consistent set of properties
	  P_ = ptmp;
	  S_ = stmp;
      }

} // namespace

#endif /* #ifndef SYNTH_INTEGRATOR*/
