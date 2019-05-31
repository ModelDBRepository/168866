/*
 *  aeif_cond_nmda_alpha.h
 *
 * Adaptive exponential integrate-and-fire model (Brette and Gerstner) with
 * alpha-shaped synaptic conductances, as implemented in aeif_cond_alpha.h 
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

#ifndef AEIF_COND_NMDA_ALPHA_H
#define AEIF_COND_NMDA_ALPHA_H

#include "config.h"

#ifdef HAVE_GSL_1_11

#include "nest.h"
#include "event.h"
#include "archiving_node.h"
#include "ring_buffer.h"
#include "connection.h"

#include "universal_data_logger.h"
#include "recordables_map.h"

//#include "analog_data_logger.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

/* BeginDocumentation
Name: aeif_cond_nmda_alpha -  Conductance based exponential integrate-and-fire
neuron model according to Brette and Gerstner (2005).

Description: aeif_cond_nmda_alpha is the adaptive exponential integrate and
fire neuron according to Brette and Gerstner (2005).  Synaptic conductances
are modelled as alpha-functions.

This implementation uses the embedded 4th order Runge-Kutta-Fehlberg solver
with adaptive stepsize to integrate the differential equation.

The membrane potential is given by the following differential equation: C
dV/dt= -g_L(V-E_L)+g_L*Delta_T*exp((V-V_T)/Delta_T)-g_e(t)(V-E_e)
-g_i(t)(V-E_i)-w +I_e

and

tau_w * dw/dt= a(V-E_L) -W

Parameters: 
The following parameters can be set in the status dictionary.

Dynamic state variables:
  V_m        double - Membrane potential in mV
  g_ex       double - Excitatory synaptic conductance in nS.
  dg_ex      double - First derivative of g_ex in nS/ms
  g_in       double - Inhibitory synaptic conductance in nS.
  dg_in      double - First derivative of g_in in nS/ms.
  w          double - Spike-adaptation current in pA.

Membrane Parameters:
  C_m        double - Capacity of the membrane in pF
  t_ref      double - Duration of refractory period in ms. 
  V_peak     double - Spike detection threshold in mV.
  V_reset    double - Reset value for V_m after a spike. In mV.
  E_L        double - Leak reversal potential in mV. 
  g_L        double - Leak conductance in nS.
  I_e        double - Constant external input current in pA.

Spike adaptation parameters:
  a          double - Subthreshold adaptation in nS.
  b          double - Spike-triggered adaptation in pA.
  Delta_T    double - Slope factor in mV
  tau_w      double - Adaptation time constant in ms
  V_t        double - Spike initiation threshold in mV (V_th can also be used for compatibility).

Synaptic parameters
  E_ex       double - Excitatory reversal potential in mV.
  tau_syn_ex double - Rise time of excitatory synaptic conductance in ms (alpha function).
  E_in       double - Inhibitory reversal potential in mV.
  tau_syn_in double - Rise time of the inhibitory synaptic conductance in ms (alpha function).

NMDA parameters
  NMDA_V_max double - Membrane voltage for the upper "shoulder" of the NMDA sigmoid 
  NMDA_V_min double - Membrane voltage for the lower "shoulder" of the NMDA sigmoid
  NMDA_gain  double - gain factor for unit sigmoid _before_ scaling 
  E_n       double - NMDA reversal potential in mV.
  tau_syn_n  double - Rise time of nmda synaptic conductance in ms (alpha function).
    
Author: Marc-Oliver Gewaltig

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

References: Brette R and Gerstner W (2005) Adaptive Exponential Integrate-and-Fire Model as 
            an Effective Description of Neuronal Activity. J
            Neurophysiol 94:3637-3642

SeeAlso: iaf_cond_alpha
*/

namespace nest
{
  class Network;

  class aeif_cond_nmda_alpha:
  public Archiving_Node
  {
    
  public:        
    
    aeif_cond_nmda_alpha();
    aeif_cond_nmda_alpha(const aeif_cond_nmda_alpha&);
    ~aeif_cond_nmda_alpha();

    /**
     * Import sets of overloaded virtual functions.
     * We need to explicitly include sets of overloaded
     * virtual functions into the current scope.
     * According to the SUN C++ FAQ, this is the correct
     * way of doing things, although all other compilers
     * happily live without.
     */
/* #ifndef IS_BLUEGENE
   using Node::check_connection;
#endif */
    using Node::connect_sender;
    using Node::handle;

    port check_connection(Connection&, port);
    
    void handle(SpikeEvent &);
    void handle(CurrentEvent &);
    void handle(DataLoggingRequest &); 
    
    port connect_sender(SpikeEvent &, port);
    port connect_sender(CurrentEvent &, port);
    port connect_sender(DataLoggingRequest &, port);

    void get_status(DictionaryDatum &) const;
    void set_status(const DictionaryDatum &);
    
  private:
    
    void init_node_(const Node& proto);
    void init_state_(const Node& proto);
    void init_buffers_();
    void calibrate();
    void update(Time const &, const long_t, const long_t);
	    friend class RecordablesMap<aeif_cond_nmda_alpha>;
	    friend class UniversalDataLogger<aeif_cond_nmda_alpha>;
  public:

    // ---------------------------------------------------------------- 

    /** 
     * Independent parameters of the model. 
     * These parameters must be passed to the iteration function that
     * is passed to the GSL ODE solvers. Since the iteration function
     * is a C++ function with C linkage, the parameters can be stored
     * in a C++ struct with member functions, as long as we just pass
     * it by void* from C++ to C++ function. The struct must be public,
     * though, since the iteration function is a function with C-linkage,
     * whence it cannot be a member function of iaf_cond_nmda_alpha.
     * @note One could achieve proper encapsulation by an extra level
     *       of indirection: Define the iteration function as a member
     *       function, plus an additional wrapper function with C linkage.
     *       Then pass a struct containing a pointer to the node and a
     *       pointer-to-member-function to the iteration function as void*
     *       to the wrapper function. The wrapper function can then invoke
     *       the iteration function on the node (Stroustrup, p 418). But
     *       this appears to involved, and the extra indirections cost.
     */
    struct Parameters_ {
      double_t V_peak_;     //!< Spike detection threshold in mV
      double_t V_reset_;    //!< Reset Potential in mV
      double_t t_ref_;      //!< Refractory period in ms

      double_t g_L;         //!< Leak Conductance in nS
      double_t C_m;         //!< Membrane Capacitance in pF
      double_t E_ex;        //!< Excitatory reversal Potential in mV
      double_t E_in;        //!< Inhibitory reversal Potential in mV
      double_t E_L;         //!< Leak reversal Potential (aka resting potential) in mV
      double_t Delta_T;     //!< Slope faktor in ms.
      double_t tau_w;       //!< adaptation time-constant in ms.
      double_t a;           //!< Subthreshold adaptation in nS.
      double_t b;           //!< Spike-triggered adaptation in pA
      double_t V_th;        //!< Spike threshold in mV.
      double_t t_ref;       //!< Refractory period in ms.
      double_t tau_syn_ex;  //!< Excitatory synaptic rise time.
      double_t tau_syn_in;  //!< Excitatory synaptic rise time.
      double_t I_e;         //!< Intrinsic current in pA.

      double_t N_V_max;      //!< NMDA max membrane voltage 
      double_t N_V_min;	  //!< NMDA min membrane voltage
      double_t N_gain;	  //!< NMDA sigmoid gain
      double_t E_n;        //!< NMDA reversal Potential in mV
      double_t tau_syn_n;   //!< NMDA synaptic rise time.
  
      /** 
       * External input current from CurrentEvents.
       * This is not a parameter but a variable. It is still placed here, since
       * it needs to be passed to the iteration function. We thus avoid the need
       * of an additional wrapper structure. It is not revealed or manipulateable.
       */
      double_t I_stim;      //!< External Stimulus in pA
  
      Parameters_();  //!< Sets default parameter values

      void get(DictionaryDatum&) const;  //!< Store current values in dictionary
      void set(const DictionaryDatum&);  //!< Set values from dicitonary
    };
  
  /**
   * Enumeration identifying elements in state array State_::y_.
   * The state vector must be passed to GSL as a C array. This enum
   * identifies the elements of the vector. It must be public to be
   * accessible from the iteration function.
   */  
  enum Statevars
  {
    V_M   = 0,
    DG_EX    ,  // 1
    G_EX     ,  // 2
    DG_IN    ,  // 3
    G_IN     ,  // 4
    W        ,  // 5
    DG_N    ,   // 6
    G_N     ,   // 7
    NSTATES
  };
    
    
  private:
    // ---------------------------------------------------------------- 

    /**
     * State variables of the model.
     * @note Copy constructor and assignment operator required because
     *       of C-style array.
     */
    struct State_ {
      double_t y_[NSTATES];  //!< neuron state, must be C-array for GSL solver
      int_t    r_;           //!< number of refractory steps remaining

      State_(const Parameters_&);  //!< Default initialization
      State_(const State_&);
      State_& operator=(const State_&);

      void get(DictionaryDatum&) const;
      void set(const DictionaryDatum&, const Parameters_&);
    };    

    // ---------------------------------------------------------------- 

    /**
     * Buffers of the model.
     */
    struct Buffers_ {
	Buffers_(aeif_cond_nmda_alpha&);
	Buffers_(const Buffers_&, aeif_cond_nmda_alpha&);
//	Buffers_(); //!<Sets buffer pointers to 0
      /** buffers and sums up incoming spikes/currents */
      RingBuffer spike_exc_;
      RingBuffer spike_inh_;
      RingBuffer spike_nmda_;
      RingBuffer currents_;
	UniversalDataLogger<aeif_cond_nmda_alpha> logger_;
      
//	AnalogDataLogger<PotentialRequest>           potentials_;
//      AnalogDataLogger<SynapticConductanceRequest> conductances_;
//      AnalogDataLogger<AeifWRequest>               aeif_ws_;

      /** GSL ODE stuff */
      gsl_odeiv_step*    s_;    //!< stepping function
      gsl_odeiv_control* c_;    //!< adaptive stepsize control function
      gsl_odeiv_evolve*  e_;    //!< evolution function
      gsl_odeiv_system   sys_;  //!< struct describing system
      
      // IntergrationStep_ should be reset with the neuron on ResetNetwork,
      // but remain unchanged during calibration. Since it is initialized with
      // step_, and the resolution cannot change after nodes have been created,
      // it is safe to place both here.
      double_t step_;           //!< step size in ms
      double   IntegrationStep_;//!< current integration time step, updated by GSL
    };

     // ---------------------------------------------------------------- 

     /**
      * Internal variables of the model.
      */
     struct Variables_ { 
      /** initial value to normalise excitatory synaptic conductance */
      double_t g0_ex_; 
   
      /** initial value to normalise inhibitory synaptic conductance */
      double_t g0_in_;    

      /** initial value to normalise NMDA synaptic conductance */
      double_t g0_n_; 

      int_t    RefractoryCounts_;
     };

    // Access functions for UniversalDataLogger -------------------------------
    
    //! Read out state vector elements, used by UniversalDataLogger
    template <Statevars elem>
    double_t get_y_elem_() const { return S_.y_[elem]; }
    // ---------------------------------------------------------------- 

    Parameters_ P_;
    State_      S_;
    Variables_  V_;
    Buffers_    B_;
    static RecordablesMap<aeif_cond_nmda_alpha> recordablesMap_;
  };

  inline  
  port aeif_cond_nmda_alpha::check_connection(Connection& c, port receptor_type)
  {

    //  printf ("\tCheckSynapse: %d\n", receptor_type);
    SpikeEvent e;
    e.set_sender(*this);
    c.check_event(e);
    return c.get_target()->connect_sender(e, receptor_type);
  }

  inline
  port aeif_cond_nmda_alpha::connect_sender(SpikeEvent&, port receptor_type)
  {
    // 0: normal synapse, 1: NMDA synapse
    //  printf ("\n\tSynapse: %d \n\n", receptor_type);
    if (receptor_type >= 2)
      throw UnknownReceptorType(receptor_type, get_name());
    return receptor_type;
  }
 
  inline
  port aeif_cond_nmda_alpha::connect_sender(CurrentEvent&, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    return 0;
  }
 
  inline
  port aeif_cond_nmda_alpha::connect_sender(DataLoggingRequest& dlr, port receptor_type)
  {
    if (receptor_type != 0)
      throw UnknownReceptorType(receptor_type, get_name());
    //B_.potentials_.connect_logging_device(pr);
    //return 0;
      return B_.logger_.connect_logging_device(dlr, recordablesMap_);
  }


  inline
  void aeif_cond_nmda_alpha::get_status(DictionaryDatum &d) const
  {
    P_.get(d);
    S_.get(d);
    Archiving_Node::get_status(d);

    (*d)[names::recordables] = recordablesMap_.get_list();
  }

  inline
  void aeif_cond_nmda_alpha::set_status(const DictionaryDatum &d)
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

#endif //HAVE_GSL
#endif //AEIF_COND_NMDA_ALPHA_H
