"""Component models for the low temperature district heating network."""
# Original Components from COMANDO Examples: Dominik Hering, Marco Langiu
# https://jugit.fz-juelich.de/iek-10/public/optimization/comando
# Modified by: Mauricio Carcamo Medel
from re import A
from comando.core import Component, BINARY, System, INTEGER
import numpy as np
from sympy import *
# from tqdm import tgrange

#######################################################
#####################   SOURCES   #####################
#######################################################

class DummySink(Component):
    """A dummy sink that serves as a sink for a massflow
    """

    def __init__(self, label):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################

        #######################################################################
        # Operational Variables
        #######################################################################
        use = self.make_operational_variable('use', bounds=(0, None),
                                             init_val=0)
        self.add_input('IN', use)

class DummySource(Component):
    """A dummy resource that serves as a source for an arbitrary commodity.
    """

    def __init__(self, label):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        price = self.make_parameter('price',0)

        #######################################################################
        # Operational Variables
        #######################################################################
        use = self.make_operational_variable('use', bounds=(0, None),
                                             init_val=0)
        self.add_output('OUT', use)

        self.add_expression('variable_costs', price * use)

class FreeMassflowSource(Component):
    """A resource that serves as a source for waste heat in T_flow and mass flow cond.
    """
    
    maxQ = 20 #[kW]
    max_mflow = 10000
    max_Pel = 25000
    def __init__(self, label):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        # Temperatures
        Tmax = self.make_parameter('T_max', 273.15 + 100)
        Tmin = self.make_parameter('T_min', 273.15 + 0)
        tflow = self.make_design_variable('tflow', bounds=(273.15,273.15+100))
        Tset = self.make_parameter('Tset', 273.15+5)
        self.add_eq_constraint(tflow,Tset, 'tflow_set')
        T_price = self.make_parameter('T_price', 0.01)

        # Mass flow rate
        mflow = self.make_operational_variable('mflow', bounds=(0, None),init_val=0)
        self.add_output(identifier='OUT_mflow', expr=mflow)
        self.add_expression('variable_costs', mflow * T_price)



#######################################################
#####################  NETWORK   ######################
#######################################################


class NetworkPipe(Component):
    """A network component. Collects and distributes energy flows.

    """

    pump_max = 4000 # [kW]
    max_mflow = 4500 # [kg/s]
    def __init__(self, label, type = 'simple',n_seg = 10):

        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        #pipe parameters
        rho_pipe = self.make_parameter('rho_pipe', 8000) # [kg/m^3] - steel - find for PE
        cp_pipe = self.make_parameter('cp_pipe', 0.5) # [kJ/Kg/K] - steel - find for PE
        #water parameters
        cp = self.make_parameter('cp', 4.184) # [kJ/Kg/K] - water = 20° C
        rho = self.make_parameter('rho', 998) # [kg/m^3] - water = 20° C
        # Ground temperature
        t_gr = self.make_parameter('T_g',273.15+10)  # [K]

        # pipe data
        pipe_cost = self.make_parameter('pipe_cost',648) # $ per m for DN400
        length = self.make_parameter('length',100) # [m]
        D_i = self.make_parameter('D_i',0.400) # [m]
        fric =  self.make_parameter('friction_factor',0.01229) #[-]
        U = self.make_parameter('U', 0.4514) # [W/m/k]
        thi = self.make_parameter('thi', 0.0063) # [m]
        m_max = self.make_parameter('m_max',389.629) #[kg/s]
        v_max = self.make_parameter('v_max',3.205391) #[m/s]
        
        # Pump efficiency
        eta_pump = self.make_parameter('eta_pump',0.72)

        #######################################################################
        # Design Variables
        #######################################################################
        Area_w = D_i**2*np.pi/4 # pipe inner area [m^2]
        Vol_w = Area_w*length #Water volume
        Vol_s = ((D_i+2*thi)**2-(D_i)**2)/4*np.pi*length #steel volume
        

        #######################################################################
        # Operational Variables
        #######################################################################
        m_flow = self.make_operational_variable("m_flow", bounds=(0, self.pump_max))
        b_op = self.make_operational_variable("b_op", domain=BINARY)
        self.add_le_constraint(m_flow, self.max_mflow * b_op, 'bop_mflow')
        self.add_input(identifier='IN_mflow', expr=m_flow)
        self.add_output(identifier='OUT_mflow', expr=m_flow)

        t_init = self.make_design_variable('t_init',bounds=(273.15,273.15+100))

        # Input temperature
        
        t_in = self.make_operational_variable('t_in',
                                                bounds=(273.15, 373.15))

        # Temperature with losses
        t_out= self.make_operational_variable('t_out',
                                                 bounds=(273.15, 373.15))


        if type == 'simple':
            self.add_eq_constraint(t_out,t_in,'constant_t')
        
        elif type == 'simple_losses':
            self.add_eq_constraint(m_max*0.453*cp * 1000 * (t_in - t_out), 
                    (t_in - t_gr) * length * U,
                    'pipe_losses')

        elif type == 'pipe_losses':
            a = self.make_parameter('a_hv',0)
            b = self.make_parameter('b_hv',0) #[kg/s]
            p_m = self.make_parameter('p_m',0)
            self.add_le_constraint(m_flow, m_max , 'massflow_tracker_le')

            self.add_ge_constraint(m_flow, p_m, 'min_mflow')
            ## defining auxiliary variable
            w_bl_sup = self.max_mflow*373.15 # maximum value of bilinear product
            w_bl_low = 0 # minimum value (no water flow)
            w_bl = self.make_operational_variable('w_bl', bounds =(w_bl_low, w_bl_sup))
            self.add_eq_constraint(w_bl, m_flow*(t_in), 'bilinear multiplication')           
            self.add_eq_constraint(t_out,
                    (w_bl*a - m_flow*a*t_gr + b*(t_in - t_gr) + t_gr),
                    'pipe_losses')


        elif type == 'extensive':
            thermal_cap = self.add_expression('thermal_cap',(rho*cp)*Area_w+(rho_pipe*cp_pipe)*Area_w*Vol_s/Vol_w)
            thermal_delay = self.add_expression('thermal_delay',cp*m_flow*(t_out-t_in)/length)

            thermal_loss = self.add_expression('thermal_loss',(U/1000*length)*Area_w/Vol_w*(t_out-t_gr)) # U from W/mK to kW/mK
            t_out_change = (thermal_delay*-1-thermal_loss)/thermal_cap*3600 # test
            self.declare_state(t_out, t_out_change, t_init)
        else:
            raise Exception("type of the pipe must be simple or extensive")

        #circulation pump

        p_pump = self.make_operational_variable('P_pump', bounds=(0, self.pump_max))
        self.add_le_constraint(p_pump, self.pump_max * b_op, 'P_HP_op_c')

        for i in range(n_seg): #linearization of massflow - power dependance params a and b must be updated accordingly after model instance generation
            a = 'a'+str(i)
            b = 'b'+str(i)
            exec('%s = self.make_parameter("%s",%f)'%(a,a,0))
            exec('%s = self.make_parameter("%s",%f)'%(b,b,0))
            exec('self.add_le_constraint((m_flow*%s+%s)*length/1000,p_pump,f"power_pump_%s_%s")' %(a,b,str(i),str(n_seg)))


        # circulation pump connected to the grid
        self.add_input(identifier='IN_P_pump', expr=p_pump)

        #######################################################################
        # Investment Costs
        #######################################################################
        inv_costs = length*pipe_cost/2 # Pipe length [m] * cost per metre [$/m] / 2 (1 of the 2 pipes)
        self.add_expression('investment_costs', inv_costs)

#######################################################
###################  TRANSFORMERS   ###################
#######################################################

class HeatSourceDecentral(Component):
    """A model for a generic decentral heat source.

    Using heat connector with efficiency of the component as a parameter

    Parameters
    ----------

    efficiency:
        Efficiency of component - unit of heat per unit of fuel(gas/elec)
    n:
        number of parallel units
    p_spec:
        Specific price of component  [€/kW]
    p_fix:
        Fix price of component [€]

    Design Variables
    ----------------
    b_build:
        Boolean for build decision of HS

    Connectors
    ----------
    OUT_mflow:
        Mass flow rate
    IN_P:
        Input for power consumption

    Expressions
    -----------
    investment_costs:    Investment cost of Heat source
    """

    Q_max = 20  # [kW] - may be modified

    def __init__(self, label):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        spec_price = self.make_parameter('p_spec', 100)  # €/kW
        fix_price = self.make_parameter('p_fix', 100)  # €
        main_price = self.make_parameter('p_main',1)        
        eff = self.make_parameter('efficiency', value=1)
        n = self.make_parameter('n', value=1)

        qdot_design = self.make_design_variable('Qdot_max',
                                                bounds=(0, self.Q_max))
        b_build = self.make_design_variable('b_build', domain=BINARY)
        self.add_le_constraint(qdot_design, b_build * self.Q_max, 'Qdot_max')
        #######################################################################
        # Operational Variables
        #######################################################################
        q_dot = self.make_operational_variable('q_dot',
                                                bounds=(0, self.Q_max))
        self.add_le_constraint(q_dot, b_build * qdot_design, 'b_qdot')
        self.add_output('OUT_qdot', q_dot * n)
        p_in = q_dot
        self.add_le_constraint(p_in, qdot_design, 'p_max')
        self.add_input('IN_P', p_in * n / eff)

        # Investment costs
        inv_costs = (spec_price * qdot_design + fix_price) * n * b_build
        self.add_expression('investment_costs', inv_costs)
        main_costs = (main_price) * n * b_build
        self.add_expression('fixed_costs', main_costs)

class ASHP(Component):
    """Quadratic Heat pump model for ASHP/RCAC. Reversible pump is considered.
    """

    Q_max = 20  # [kW]
    mdot_max = 100  # [kg/s]

    def __init__(self, label,
                 t_out_b_h,t_out_b_c, mode='reversible'):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        eta_nom_h = self.make_parameter('eta_h', 0.4) #heating efficiency
        eta_nom_c = self.make_parameter('eta_c', 0.4) # cooling efficiency
        cp = self.make_parameter('cp', 4.184)
        n = self.make_parameter('n', 1)
        T_amb = self.make_parameter('T_amb',10+273.15)
        dT_a_h = self.make_parameter('dT_a_h',2.5)
        dT_a_c = self.make_parameter('dT_a_c',-2.5)
        spec_price = self.make_parameter('p_spec', 100)  # €/kW
        fix_price = self.make_parameter('p_fix', 100)  # €
        main_price = self.make_parameter('p_main',1) # €
        max_cop_h = self.make_parameter('max_cop_h', 5)
        max_cop_c = self.make_parameter('max_cop_c', 5)
        # dT a
        t_in_a = T_amb  # cold temperature [K]
        t_out_a_h = T_amb - dT_a_h # delta in evaporator - heating mode
        t_out_a_c = T_amb - dT_a_c # delta in condenser - cooling mode
        #######################################################################
        # Design Variables
        #######################################################################
        qdot_design = self.make_design_variable('Qdot_design',
                                                bounds=(0, self.Q_max))
        b_build = self.make_design_variable('b_build', domain=BINARY,
                                            bounds=(0, 1))
        self.add_le_constraint(qdot_design, self.Q_max * b_build, 'HP_Q_max')

        #######################################################################
        # Operational Variable
        #######################################################################


        b_h = self.make_operational_variable('b_h', domain=BINARY)
        b_c = self.make_operational_variable('b_c',domain=BINARY)

        if mode == 'reversible':
            self.add_le_constraint(b_c, 1 - b_h,'Heating_or_cooling')
        elif mode == 'heat_only':
            self.add_le_constraint(b_h, 1, 'heating_mode')
            self.add_le_constraint(b_c, 0, 'no_cooling')
        elif mode == 'cool_only':
            self.add_le_constraint(b_c, 1, 'cooling_mode')
            self.add_le_constraint(b_h, 0, 'no_heating')
        else:
            raise TypeError("specify either reversible, heat_only or cool_only for ashp operation")

        self.add_le_constraint(b_h, b_build, 'HP_heating_operation')
        self.add_le_constraint(b_c, b_build, 'HP_cooling_operation')

        ######## heating ########

        # Evaporator (side a)
        mflow_a_h = self.make_operational_variable('mflow_a_h',
                                                    bounds=(0, self.mdot_max))
        self.add_le_constraint(mflow_a_h, b_h * self.mdot_max,
                            'HP_mflow_a_h')
        self.add_le_constraint(mflow_a_h, b_build * self.mdot_max,
                            'HP_mflow_a_h_build')
        self.add_input(identifier='IN_mflow_a_h', expr=mflow_a_h * n)

        qdot_a_h = mflow_a_h * (t_in_a - t_out_a_h) * cp

        # Condenser (side b)

        qdot_b_h = self.make_operational_variable('qdot_b_h',
                                                    bounds=(0, self.Q_max))
        self.add_le_constraint(qdot_b_h, b_h * qdot_design,
                            'HP_qdot_b_h_op')
        self.add_le_constraint(qdot_b_h, b_build * qdot_design,
                            'HP_qdot_b_h_build')
        self.add_output(identifier='OUT_qdot_h', expr=qdot_b_h * n)

        # HP - heating
        COP_h = self.make_operational_variable('COP_h', bounds=(0, 10))
        self.add_le_constraint(COP_h,max_cop_h *b_h, 'Maximum_COP_heating')  
        self.add_le_constraint(qdot_b_h, self.Q_max * b_h, 'Qdot_HP_op_h')

        p_hp_h = self.make_operational_variable('P_HP_h', bounds=(0, self.Q_max))
        self.add_le_constraint(p_hp_h, self.Q_max * b_h, 'P_HP_op_h')

        # n_HP Heat pumps are connected to the grid in parallel
        self.add_input(identifier='IN_P_h', expr=p_hp_h * n)
        self.add_le_constraint(qdot_b_h, qdot_design, 'output_limit_h')

        
        
        self.add_eq_constraint(p_hp_h * COP_h,
                            qdot_b_h,
                            'input_output_relation_heating')
        self.add_le_constraint(COP_h*(t_out_b_h - t_out_a_h), t_out_b_h * eta_nom_h*b_h, 'COP_h_calc')
        self.add_le_constraint(COP_h,max_cop_h, 'COP_h_max')
        self.add_eq_constraint(qdot_a_h + p_hp_h, qdot_b_h, 'energy_balance_h')

        ######## cooling ######## 

        # Condenser (side a)
        mflow_a_c = self.make_operational_variable('mflow_a_c',
                                                    bounds=(0, self.mdot_max))
        self.add_le_constraint(mflow_a_c, b_c * self.mdot_max,
                            'HP_mflow_a_c')
        self.add_le_constraint(mflow_a_c, b_build * self.mdot_max,
                            'HP_mflow_a_c_build')
        self.add_input(identifier='IN_mflow_a_c', expr=mflow_a_c * n)

        qdot_a_c = mflow_a_c * -(t_in_a - t_out_a_c) * cp ## injected energy in outer fluid

        # Evaporator (side b)

        qdot_b_c = self.make_operational_variable('qdot_b_c',
                                                    bounds=(0, self.Q_max))
        self.add_le_constraint(qdot_b_c, b_c * qdot_design,
                            'HP_qdot_b_c_op')
        self.add_le_constraint(qdot_b_c, b_build * qdot_design,
                            'HP_qdot_b_c_build')
        self.add_output(identifier='OUT_qdot_c', expr=qdot_b_c * n)

        # HP - cooling
        COP_c = self.make_operational_variable('COP_c', bounds=(0, 10))
        self.add_le_constraint(qdot_b_c, self.Q_max * b_c, 'Qdot_HP_op_c')

        p_hp_c = self.make_operational_variable('P_HP_c', bounds=(0, self.Q_max))
        self.add_le_constraint(p_hp_c, self.Q_max * b_c, 'P_HP_op_c')
                
        # n_HP Heat pumps are connected to the grid in parallel
        self.add_input(identifier='IN_P_c', expr=p_hp_c * n)
        self.add_le_constraint(qdot_b_c, qdot_design, 'output_limit_c')

        self.add_eq_constraint(p_hp_c * COP_c,
                            qdot_b_c, 
                            'input_output_relation_cooling')
        self.add_le_constraint(COP_c *  (t_out_a_c - t_out_b_c),
                        t_out_b_c * eta_nom_c, 
                            'COP_c_calc')                                            
        self.add_le_constraint(COP_c,max_cop_c, 'COP_c_max')
        self.add_eq_constraint(qdot_b_c + p_hp_c,qdot_a_c, 'energy_balance_c')        

        #######################################################################
        # Misc
        #######################################################################
        # Investment costs
        inv_costs = (spec_price * qdot_design + b_build * fix_price) * n
        self.add_expression('investment_costs', inv_costs)
        main_costs = (main_price) * n * b_build
        self.add_expression('fixed_costs', main_costs)

class DecenHeatPump(Component):
    """Quadratic Heat pump model.
    """

    Q_max = 20  # [kW]
    mdot_max = 10  # [kg/s]

    def __init__(self, label,
                 t_out_b_h,t_out_b_c,b_connect = None, dt_max = 15, mode = 'reversible'):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        eta_nom_h = self.make_parameter('eta_h', 0.55) #heating efficiency
        eta_nom_c = self.make_parameter('eta_c', 0.50) # cooling efficiency
        cp = self.make_parameter('cp', 4.184)
        n = self.make_parameter('n', 1)
        pinch_dt = self.make_parameter('pinch_dt',0.0) ### Internal HP temperatures
        spec_price = self.make_parameter('p_spec', 100)  # €/kW
        fix_price = self.make_parameter('p_fix', 100)  # €
        main_price = self.make_parameter('p_main',1) # €
        max_cop_h = self.make_parameter('max_cop_h', 7)
        max_cop_c = self.make_parameter('max_cop_c', 6)
        annuity_factor = self.make_parameter('ann_factor',1)
        #######################################################################
        # Design Variables
        #######################################################################
        qdot_design = self.make_design_variable('Qdot_design',
                                                bounds=(0, self.Q_max))
        b_build = self.make_design_variable('b_build', domain=BINARY,
                                            bounds=(0, 1))
        self.add_le_constraint(qdot_design, self.Q_max * b_build, 'HP_Q_max')

        if b_connect is not None:
            self.add_le_constraint(b_build, b_connect, 'connection')


        #######################################################################
        # Operational Variable
        #######################################################################
        dT_a_h = self.make_operational_variable('dT_a_h',bounds=(0, dt_max))
        dT_a_c = self.make_operational_variable('dT_a_c',bounds=(-1 * 25, 0)) ## feasible in all cases for the given conditions
        b_h = self.make_operational_variable('b_h', domain=BINARY)
        b_c = self.make_operational_variable('b_c',domain=BINARY)
        if mode == 'reversible':
            self.add_le_constraint(b_c, 1 - b_h,'Heating_or_cooling')
        elif mode == 'heat_only':
            self.add_le_constraint(b_h, 1, 'heating_mode')
            self.add_le_constraint(b_c, 0, 'no_cooling')
        elif mode == 'cool_only':
            self.add_le_constraint(b_c, 1, 'cooling_mode')
            self.add_le_constraint(b_h, 0, 'no_heating')
        else:
            raise TypeError("specify either reversible, heat_only or cool_only for ashp operation")
        
        self.add_le_constraint(b_h, b_build, 'HP_heating_operation')
        self.add_le_constraint(b_c, b_build, 'HP_cooling_operation')

        
        ######## temperature in side A #########
        # dT a
        t_in_a = self.make_operational_variable('t_in_a', bounds=(273.15,273.15+100))
        t_out_a_h = self.make_operational_variable('t_out_a_h', bounds=(273.15,273.15+100))
        t_out_a_c = self.make_operational_variable('t_out_a_c', bounds=(273.15,273.15+100))

        self.add_eq_constraint(t_out_a_h, t_in_a - dT_a_h, 'side_a_dT_heating')
        self.add_eq_constraint(t_out_a_c, t_in_a - dT_a_c, 'side_a_dT_cooling')

    

        ######## heating ########



        # Evaporator (side a)
        mflow_a_h = self.make_operational_variable('mflow_a_h',
                                                    bounds=(0, self.mdot_max))
        self.add_le_constraint(mflow_a_h, b_h * self.mdot_max,
                            'HP_mflow_a_h_op')
        self.add_le_constraint(mflow_a_h, b_build * self.mdot_max,
                            'HP_mflow_a_h_des')
        self.add_input(identifier='IN_mflow_a_h', expr=mflow_a_h * n)
        
        self.add_output(identifier='OUT_mflow_a_h', expr=mflow_a_h * n)


        qdot_a_h = self.make_operational_variable('qdot_a_h',
                                                    bounds=(0, self.Q_max))
        self.add_le_constraint(qdot_a_h, b_h * qdot_design,
                            'HP_qdot_a_h_op')
        self.add_le_constraint(qdot_a_h, b_build * qdot_design,
                            'HP_qdot_a_h_des')
        self.add_eq_constraint(qdot_a_h,mflow_a_h * (dT_a_h) * cp, 'heating_side_a_q')

        # Condenser (side b)
        #qdot
        qdot_b_h = self.make_operational_variable('qdot_b_h',
                                                    bounds=(0, self.Q_max))
        self.add_le_constraint(qdot_b_h, b_h * qdot_design,
                            'HP_qdot_b_h_op')
        self.add_le_constraint(qdot_b_h, b_build * qdot_design,
                            'HP_qdot_b_h_des')
        self.add_output(identifier='OUT_qdot_h', expr=qdot_b_h * n)


        # HP - heating
        COP_h = self.make_operational_variable('COP_h', bounds=(0, 10))
        self.add_le_constraint(COP_h,max_cop_h *b_h, 'Maximum_COP_heating')  
        self.add_le_constraint(qdot_b_h, self.Q_max * b_h, 'Qdot_HP_op_h')

        p_hp_h = self.make_operational_variable('P_HP_h', bounds=(0, self.Q_max))
        self.add_le_constraint(p_hp_h, self.Q_max * b_h, 'P_HP_op_h')
        t_cond_h = t_out_b_h + pinch_dt
        t_evap_h = t_out_a_h - pinch_dt
        self.add_eq_constraint(p_hp_h * COP_h,
                            qdot_b_h,
                            'input_output_relation_heating')
        self.add_le_constraint(COP_h*(t_cond_h - t_evap_h), t_cond_h * eta_nom_h*b_h, 'COP_h_calc')
        self.add_le_constraint(COP_h,max_cop_h, 'COP_h_max')
        self.add_eq_constraint(qdot_a_h + p_hp_h, qdot_b_h, 'energy_balance_h')


        # n_HP Heat pumps are connected to the grid in parallel
        self.add_input(identifier='IN_P_h', expr=p_hp_h * n)
        self.add_le_constraint(qdot_b_h, qdot_design, 'output_limit_h')


        ######## cooling ########

        # Condenser (side a)
        mflow_a_c = self.make_operational_variable('mflow_a_c',
                                                    bounds=(0, self.mdot_max))
        self.add_le_constraint(mflow_a_c, b_c * self.mdot_max,
                            'HP_mflow_a_c_op')
        self.add_le_constraint(mflow_a_c, b_build * self.mdot_max,
                            'HP_mflow_a_c_des')
        self.add_input(identifier='IN_mflow_a_c', expr=mflow_a_c * n)
        self.add_output(identifier='OUT_mflow_a_c', expr=mflow_a_c * n)

        qdot_a_c = self.make_operational_variable('qdot_a_c',
                                                    bounds=(0,self.Q_max))
        self.add_le_constraint(qdot_a_c, self.Q_max * b_c, 'Qdot_HP_op_c')
        self.add_le_constraint(qdot_a_c, self.Q_max * b_build, 'Qdot_HP_des_c')
        self.add_eq_constraint(qdot_a_c,mflow_a_c * -1*(dT_a_c) * cp, 'cooling_side_a_c') # injected energy in outer fluid

        # Evaporator (side b)
        #qdot
        qdot_b_c = self.make_operational_variable('qdot_b_c',
                                                    bounds=(0, self.Q_max))
        self.add_le_constraint(qdot_b_c, b_c * qdot_design,
                            'HP_qdot_b_c_op')
        self.add_le_constraint(qdot_b_c, b_build * qdot_design,
                            'HP_qdot_b_c_des')
        self.add_output(identifier='OUT_qdot_c', expr=qdot_b_c * n)

        # HP - cooling
        
        COP_c = self.make_operational_variable('COP_c', bounds=(0, 10))
        self.add_le_constraint(qdot_b_c, self.Q_max * b_c, 'Qdot_HP_op_c')
        self.add_le_constraint(qdot_b_c, self.Q_max * b_c, 'Qdot_HP_op_c')
        self.add_le_constraint(qdot_b_c, self.Q_max * b_build, 'Qdot_HP_des_c')

        p_hp_c = self.make_operational_variable('P_HP_c', bounds=(0, self.Q_max))
        self.add_le_constraint(p_hp_c, self.Q_max * b_c, 'P_HP_op_c')
        t_cond_c = t_out_a_c + pinch_dt
        t_evap_c = t_out_b_c - pinch_dt
        t_op = self.make_operational_variable('temp_op_c', domain=BINARY)
        self.add_le_constraint(t_out_a_c, (t_out_b_c) + 100 * t_op, 'cooling_input_temperature_1')
        self.add_ge_constraint(t_out_a_c, t_out_b_c * t_op, 'cooling_input_temperature_2')

        self.add_eq_constraint(p_hp_c * COP_c,
                            qdot_b_c, 
                            'input_output_relation_cooling')
        self.add_le_constraint(COP_c *  (t_cond_c - t_evap_c),
                t_evap_c * eta_nom_c, 
                    'COP_c_calc')                                  
        self.add_le_constraint(COP_c,max_cop_c, 'COP_c_max')
            
        self.add_eq_constraint(qdot_b_c + p_hp_c,qdot_a_c, 'energy_balance_c') 

        # n_HP Heat pumps are connected to the grid in parallel
        self.add_input(identifier='IN_P_c', expr=p_hp_c * n)
        self.add_le_constraint(qdot_b_c, qdot_design, 'output_limit_c')
    

        #######################################################################
        # Misc
        #######################################################################
        # Investment costs
        inv_costs = (spec_price * qdot_design +  fix_price) * n * b_build
        self.add_expression('investment_costs', inv_costs)
        main_costs = (spec_price/annuity_factor * qdot_design +  fix_price/annuity_factor) * n * b_build * 0.025
        self.add_expression('fixed_costs', main_costs)

class DecenHeatExchanger(Component):
    """Heat exchanger model.

    Based on the following temperatures:
    - t_in_a:
        Incoming temperature at side a [K]
    - t_out_a:
        Outgoing temperature at side a [K]
    - t_out_b:
        Outgoing temperature at side b [K]
    - t_in_b:
        Incoming temperature at side b [K]

    Parameters
    ----------
    n:
        Integer, Number of HP instances, default = 1
    cp:
        Heat capacity of fluid, default = 4.184 kJ/(kgK)

    Design Variables
    ----------------
    b_build:
        Boolean for build decision of HX
    Qdot_design:
        Maximum heating capacity [kW]

    Operational Variables
    ---------------------
    b_op:
        Boolean operational variable.
        b_op=0 for no mass flow and no temperature constraints
        b_op=1 for heat exchanger operation
    mflow_a:
        Incoming mass flow rate at side a [kg/s]
    mflow_b:
        Outgoing mass flow rate at side b [kg/s]

    Connectors
    ----------
    IN_mflow_a:
        mass flow rate at side a [kg/s]
    OUT_mflow_b:
        mass flow rate at side b [kg/s]

    Expressions
    -----------
        investment_costs:    Investment cost of HX
    """

    Q_max = 400  # [kW]
    mdot_max = 10  # [kg/s]

    def __init__(self, label, t_out_b, b_build_hp=None):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        cp = self.make_parameter('cp', 4.184)
        n = self.make_parameter('n', 1)
        dT_a = self.make_parameter('dt_a', 8)
        pinch_dt = self.make_parameter('dt_pinch',2)
        spec_price = self.make_parameter('p_spec', 150)  # €/kW
        fix_price = self.make_parameter('p_fix', 3.75)  # €
        main_price = self.make_parameter('p_main',1) # €

        #######################################################################
        # Design Variables
        #######################################################################
        qdot_design = self.make_design_variable('Qdot_design',
                                                bounds=(0, self.Q_max))
        b_build = self.make_design_variable('b_build', domain=BINARY)
        if b_build_hp is not None:
            # Opposite of HP build decision
            self.add_le_constraint(b_build, 1 - b_build_hp)
        self.add_le_constraint(qdot_design, self.Q_max * b_build, 'Qdot_max')


        #######################################################################
        # Operational Variable
        #######################################################################
        b_op = self.make_operational_variable('b_op', domain=BINARY)
        t_in_a = self.make_operational_variable('t_in_a', bounds=(273.15,273.15+100))
        t_out_a = self.make_operational_variable('t_out_a', bounds=(273.15,273.15+100))
        # t_in_b = self.make_operational_variable('t_in_b', bounds=(273.15,273.15+100))

        self.add_le_constraint(b_op, b_build, 'operation')
        # Side a
        mflow_a = self.make_operational_variable('mflow_a',
                                                 bounds=(0, self.mdot_max))
        self.add_le_constraint(mflow_a, b_op * self.mdot_max, 'b_hx_mflow_a')
        self.add_le_constraint(mflow_a, b_build * self.mdot_max,
                               'bbuild_hx_mflow_a')
        self.add_input(identifier='IN_mflow_a', expr=mflow_a * n)
        self.add_input(identifier='IN_t_a', expr=t_in_a)
        self.add_output(identifier='OUT_t_a_h', expr=t_out_a)
        self.add_output(identifier='OUT_mflow_a', expr=mflow_a * n)
        self.add_eq_constraint(t_out_a, t_in_a - dT_a, 'side_a_dT_heating')

        qdot_in = mflow_a * (dT_a) * cp

        # Side b
        qdot_b = self.make_operational_variable('qdot_b_h',
                                                        bounds=(0, self.Q_max))
        self.add_le_constraint(qdot_b, b_op * qdot_design,
                                'HP_qdot_b_h_op')
        self.add_le_constraint(qdot_b, b_build * qdot_design,
                                'HP_qdot_b_h_des')
        self.add_output(identifier='OUT_qdot_h', expr=qdot_b * n)



        self.add_le_constraint(qdot_b, qdot_design, 'design_Limit')
        self.add_le_constraint(t_out_b * b_op, (t_in_a - pinch_dt) * b_op,
                                'max_temp_increase')

        self.add_eq_constraint(qdot_in, qdot_b, 'heat_flow')

        #######################################################################
        # Misc
        #######################################################################
        # Investment costs
        inv_costs = (spec_price * qdot_design +  fix_price) * n * b_build
        self.add_expression('investment_costs', inv_costs)
        main_costs = main_price * n * b_build
        self.add_expression('fixed_costs', main_costs)


class LinkingComponentConsumer(System):
    """Create a linking component, modeled as a COMANDO system.

    This component contains a Heat Pump (Decentralised) and a Heat Exchanger model.
    """

    def __init__(self, label,prefix,
                 t_out_b_h,t_out_b_c,dt_max):
        super().__init__(prefix+label)

        hp = DecenHeatPump(f'{prefix}WSHP_{label}',
            t_out_b_h=t_out_b_h,t_out_b_c=t_out_b_c, dt_max = dt_max)

        hx = DecenHeatExchanger(f'{prefix}HX_{label}', b_build_hp=hp['b_build'],
            t_out_b=t_out_b_h)

        comps_h = [hp, hx]
        for comp in comps_h:
            self.add(comp)



        # Side a
        self.connect('IN_mflow_a_h', [
                                    hp.IN_mflow_a_h,
                                    hx.IN_mflow_a])
        self.extend_connection('IN_mflow_a_h')        

        self.connect('OUT_mflow_a_h', [hp.OUT_mflow_a_h,
                                    hx.OUT_mflow_a])
        self.extend_connection('OUT_mflow_a_h')

        # Side b_h
        self.connect('OUT_qdot_h', [hp.OUT_qdot_h,
                                    hx.OUT_qdot_h])
        self.extend_connection('OUT_qdot_h')
        # Side b_c
        self.connect('OUT_qdot_c', [hp.OUT_qdot_c])
        self.extend_connection('OUT_qdot_c')
                                     
        self.expose_connector(hp.IN_mflow_a_c,'IN_mflow_a_c')
        self.expose_connector(hp.OUT_mflow_a_c,'OUT_mflow_a_c')
        self.expose_connector(hp.IN_P_c, 'IN_P_c')
        self.expose_connector(hp.IN_P_h, 'IN_P_h')



#######################################################
#####################   CONSUMER   ####################
#######################################################


class BESMFlowTFixDT(Component):
    """A model for a building energy system.
    """

    def __init__(self, label):
        super().__init__(label)
        #######################################################################
        # Parameters
        #######################################################################
        qdot = self.make_parameter('Qdot')  # [kW]
        self.make_parameter('T_flow')
        self.add_connector('IN_qdot', qdot)


class Consumer(System):
    """Create a consumer group, modeled as a COMANDO system.

    This component contains a linking component, two heat sources and one
    demand.
    """

    def __init__(self, label,**kwargs):
        prefix = kwargs.get('prefix','')
        dt_max = kwargs.get('dt_max',15)
        super().__init__(prefix + label)
        building_h = BESMFlowTFixDT(f'{prefix}BES_{label}_heating')
        building_c = BESMFlowTFixDT(f'{prefix}BES_{label}_cooling')


        lc = LinkingComponentConsumer(f'LC_{label}',prefix,
                                    t_out_b_h=building_h['T_flow'],t_out_b_c=building_c['T_flow'], dt_max = dt_max)


        for comp in [lc,building_h,building_c]:                                   
            self.add(comp)
    
        self.connect('qdot_to_bes_h', [
                                      lc.OUT_qdot_h,
                                      building_h.IN_qdot])

        self.connect('qdot_to_bes_c', [
                                    lc.OUT_qdot_c,
                                    building_c.IN_qdot])

        self.connect('IN_P_el', [lc.IN_P_c,
                                 lc.IN_P_h,
                                 ])

        self.expose_connector(lc.IN_mflow_a_h, 'IN_mflow_a_h')
        self.expose_connector(lc.IN_mflow_a_c, 'IN_mflow_a_c')
        self.expose_connector(lc.OUT_mflow_a_h, 'OUT_mflow_a_h')
        self.expose_connector(lc.OUT_mflow_a_c, 'OUT_mflow_a_c')