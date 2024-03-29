'''Decen case for VIC & NRW
Preprocessing, Execution and Postprocessing'''
# AUTHOR: Mauricio Carcamo Medel.

#### file - database handling and housekeeping
import os
from pathlib import Path

### Run definitions
SEED = 123
region = 'NRW' # NRW or VIC
num_typ_days = 12 ## Available: 12. None for 8760 hours analysis.
if num_typ_days == None:
    typ_day_file = None
else:
    typ_day_file = region + '_Medoid_' + str(num_typ_days) + 'd.csv'

ins = 'S3'
run = region + '_DECEN_' + ins


local_folder = 'Base_code'
## comment or un comment where appropiate

el_price_list = [] ## electricity price vector
append_base = False ## append base electricity price value to the vector
run = region + '_DECEN_' + str(num_typ_days) + 'd_elpr_010203_' + ins

Tin_vect = [x + 273.15 for x in [10,30,50]] #
len_vect = [50,100,500,1000,5000,10000] ### 
n_con = [100,1000,10000]
reprocessing = False ### if reprocessing a existing dataset, reprocessing = True, otherwise False
timelimit = False ### if reprocessing runs that reached timelimit (plus reprocessing = True)
# repro_name_list = [run] ### if empty, no reprocessing will occur
# repro_timelimit_list = [run] ### if empty, no reprocessing will occur
 
#### script import
parent_path = Path(__file__).parents[0]
os.chdir(str(parent_path))

#### Utilities functions
def find_neighbours(value, df, colname):
    exactmatch = df[df[colname] == value]
    if not exactmatch.empty:
        return exactmatch.index
    else:
        upperneighbour_ind = df[df[colname] > value][colname].idxmin()
        return upperneighbour_ind

def upper_neighbours(value, df, colname):
    upperneighbour_ind = df[df[colname] > value][colname].idxmin()
    return upperneighbour_ind

def DN_calc(dT,T,q_c,q_h,n_con,pipeset,L, cooling = 'yes'): ## calculation of transmission pipeline inner diameter from catalogue

    import numpy as np
    #heating
    cp = 4.184
    rho = 998
    eff = 0.50
    Th_hot = 55+273.15
    DN_ind = 0
    max_flow = pipeset['Max_mflow'].iloc[DN_ind]
    #Tg
    if region == 'VIC':
        Tg_h = 9.6 + 273.15
        Tg_c = 17.92 + 273.15
    elif region == 'NRW':
        Tg_h = 4.6 + 273.15
        Tg_c = 14.98+273.15
    else:
        ValueError('Region must be either NRW or VIC, if another region is being included DN_calc function must be updated accordingly')

    # insulation
    R = pipeset['Res_' + ins].iloc[DN_ind]
    MaxTout_h = (T - Tg_h) * np.exp(-L/(rho*cp*max_flow*R)) + Tg_h

    dt_vect = np.linspace(0.5,25,int(25/0.5)).tolist()
    checker_h = True
    m_max_h_vect = []
    while checker_h:
        for dt in dt_vect:
            if T-dt < 273.15:
                m_max_h_vect.append(1E8)      
            else:          
                COP_H =Th_hot/(Th_hot-(MaxTout_h-dt))*eff
                if COP_H >= 7.0:
                    COP_H = 7.0
                m_max_h_vect.append(q_h/(dt*cp)*(1-1/COP_H)*1.0*n_con)     
        if any(i<= max_flow for i in m_max_h_vect):
            checker_h = False
        else:
            DN_ind = upper_neighbours(pipeset['DN'].iloc[DN_ind],pipeset,'DN')
            max_flow = pipeset['Max_mflow'].iloc[DN_ind]   
            R = pipeset['Res_' + ins].iloc[DN_ind]        
            MaxTout_h = (T - Tg_h) * np.exp(-L/(rho*cp*max_flow*R)) + Tg_h
            m_max_h_vect = []

    #cooling
    MaxTout_c = (T - Tg_c) * np.exp(-L/(rho*cp*max_flow*R)) + Tg_c
    Tc_cold = 7 + 273.15
    
    m_max_c_vect = []
    checker_c = True
    while checker_c:
        for dt in dt_vect:    
            Tc_hot = MaxTout_c+dt
            if Tc_hot <= Tc_cold:
                COP_C = 6.0
            else:
                COP_C =Tc_cold/(Tc_hot-Tc_cold)*eff
            if COP_C >= 6.0:
                COP_C = 6.0
            m_max_c_vect.append(q_c/(dt*cp)*(1+1/COP_C)*1.0*n_con)

        if any(i<= max_flow for i in m_max_c_vect):
            checker_c = False
        else:
            DN_ind = upper_neighbours(pipeset['DN'].iloc[DN_ind],pipeset,'DN')
            max_flow = pipeset['Max_mflow'].iloc[DN_ind]   
            R = pipeset['Res_' + ins].iloc[DN_ind]        
            MaxTout_c = (T - Tg_c) * np.exp(-L/(rho*cp*max_flow*R)) + Tg_c
            m_max_c_vect = []
    return pipeset['DN'].iloc[DN_ind],MaxTout_c,max_flow

def calc_func_mult(pset): #function for model execution, folder creation and data storage

    import pickle
    import pandas as pd
    import numpy as np
    from _scripts.utilities import annuity
    # import _scripts.pyomoio as po
    from pyomo.environ import value
    
    targetdir_abs = parent_path / pset['_calc_dir'] / pset['_pset_id']
    targetdir_abs.mkdir(parents=True,exist_ok=True)

    return {
    'pset_dict': pset,
    }

def post_process_cost(pset,run,mode='pickle'): #function for post processing of the costs
    import pickle
    import pandas as pd
    import numpy as np
    from _scripts.utilities import annuity
    from pyomo.environ import value

    calc_dirs_list = [x for x in pset['_calc_dir']]
    pset_id_list = [x for x in pset['_pset_id']]
    typ_days_list = [x for x in pset['typical_days']]
    new_calc_dirs_list = [str(parent_path) + x.split(local_folder)[1] for x in pset['_calc_dir']]
    Dic_costs = {}
    for i,(typ_days,c_dir,pset_id) in enumerate(zip(typ_days_list,new_calc_dirs_list,pset_id_list)):
        targetdir = c_dir + '/' + pset_id
        if mode == 'pickle': # if the whole pickle is save for each run (base case)
            pickledir = targetdir + '/' + 'details' + '.pk'
            with open(f'{pickledir}', 'rb') as f: # details pickle
                    details = pickle.load(f)
            params = details['parameters']
            operation = details['operation']
            design = details['design']
            if 'scenarios' in details:
                scenarios = details['scenarios']
            tot_demand = pset['tot_dem'][i]
            heat_tot_demand = pset['h_dem'][i]
            cool_tot_demand = pset['c_dem'][i]

        elif mode == 'no_pickle': # if the pickle file is not save due to storage issues (sensitivity)
            params = pd.read_csv(str(targetdir) + '/' + 'parameters_inp.csv')
            operation = pd.read_csv(str(targetdir) + '/' + 'operation_var.csv')
            design = pd.read_csv(str(targetdir) + '/' + 'design_var.csv')
            heat_tot_demand = params.BES_C_cons_heating_Qdot.sum() / 1e3  # in MWh
            cool_tot_demand = params.BES_C_cons_cooling_Qdot.sum() / 1e3  # in MWh
            tot_demand = (heat_tot_demand + cool_tot_demand) # in MWh
        
        comp_list = ['HS_gas','HS_el','Pipes','WSHP_LC','Cen_HX'] ## update LC to HX eventually
        cons_n = params.filter(regex = '_cons_n').iloc[0,0]

        if typ_days == None:

            ### Operational costs - Total commodities - no typical days
            regex_price = 'el_ind_price|el_cons_price|gas_price'
            regex_commodities = 'el_ind_use|el_cons_use|gas_use'
            commodities = params.filter(regex = regex_price).join(operation.filter(regex = regex_commodities))
            total_commodities = (commodities.filter(regex='Source_el_cons_price|Source_el_cons_use').prod(axis=1)+commodities.filter(regex='Source_el_ind_price|Source_el_ind_use').prod(axis=1) + commodities.filter(like='gas').prod(axis=1)).sum()

            ## Operational costs - Op costs per components # not considering ASHP
            regex_op_costs = 'Inflow_pipe_P|Return_pipe_P|WSHP_LC_C_cons_P_HP|HS_el_C_cons_q_dot|Source_gas_use' ### using gas use instead of el_gas_Qdot because of boiler efficiency
            op_costs = params.filter(regex = regex_price).join(operation.filter(regex = regex_op_costs))
            op_costs['weights'] = 1 # all simulated days contribute in the same proportion
            op_costs_comp = pd.DataFrame(columns=comp_list)
        else:
            scenarios_dict = scenarios.to_dict()['pi']
            ### Operational costs - Total commodities
            regex_price = 'el_ind_price|el_cons_price|gas_price'
            regex_commodities = 'el_ind_use|el_cons_use|gas_use'
            commodities = params.filter(regex = regex_price).join(operation.filter(regex = regex_commodities))
            commodities.reset_index(inplace=True)
            commodities['weights'] = commodities['s'].map(scenarios_dict.get).apply(lambda x: x)
            total_commodities = (commodities.filter(regex='Source_el_cons_price|Source_el_cons_use|weights').prod(axis=1)+commodities.filter(regex='Source_el_ind_price|Source_el_ind_use|weights').prod(axis=1) + commodities.filter(like='gas').prod(axis=1)).sum()

            ## Operational costs - Op costs per components
            regex_op_costs = 'Inflow_pipe_P|Return_pipe_P|ASHP_C_cons_P_HP|WSHP_LC_C_cons_P_HP|HS_el_C_cons_q_dot|Source_gas_use' ### using gas use instead of el_gas_Qdot because of boiler efficiency
            op_costs = params.filter(regex = regex_price).join(operation.filter(regex = regex_op_costs))
            op_costs.reset_index(inplace=True)
            op_costs['weights'] = commodities['s'].map(scenarios_dict.get).apply(lambda x: x)
            op_costs_comp = pd.DataFrame(columns=comp_list)

            ## Heat losses/gains in pipes 
            # filter relevant variables (Tin-Tout,mflow and heating/cooling loads)
            losses = operation.filter(regex = 'pipe_t_in|pipe_t_out|pipe_m_flow').join(params.filter(regex='BES_C_cons_heating_Qdot|BES_C_cons_cooling_Qdot'))
            losses_col = losses.columns.tolist()# Getting variable names Inflow_pipe_m_flow 
            #obtain variable names
            variables_losses = ['Inflow_pipe_t_out','Inflow_pipe_t_in','Inflow_pipe_m_flow','Return_pipe_t_out','Return_pipe_t_in','Return_pipe_m_flow']
            v_names_loss = {key:None for key in variables_losses}
            for x in variables_losses:
                index = [i for i, s in enumerate(losses_col) if x in s]
                if len(index) == 1:
                    v_names_loss[x] = losses_col[index[0]]
                elif len(index) == 0:
                    IndexError('Variables of interest in heat gains/losses have 0 equal names')
                else:
                    IndexError('Variables of interest in heat gains/losses have repeated names (2 or more')
            #calculate losses/gains
            losses['losses_inflow'] = losses.apply(lambda row: row[v_names_loss['Inflow_pipe_m_flow']] * (row[v_names_loss['Inflow_pipe_t_out']] - row[v_names_loss['Inflow_pipe_t_in']]) * 4.184, axis = 1) # in kWh
            losses['losses_return'] = losses.apply(lambda row: row[v_names_loss['Return_pipe_m_flow']] * (row[v_names_loss['Return_pipe_t_out']] - row[v_names_loss['Return_pipe_t_in']]) * 4.184, axis = 1) # in kWh
            losses['tot_losses'] = losses.apply(lambda row: row['losses_inflow'] + row['losses_return'], axis = 1) # in kWh
            losses['tot_generated'] = losses.apply(lambda row: row[v_names_loss['Return_pipe_m_flow']] * (row[v_names_loss['Inflow_pipe_t_out']] - row[v_names_loss['Return_pipe_t_in']]) * 4.184, axis = 1) # in kWh
            losses['tot_generated_source'] = losses.apply(lambda row: row[v_names_loss['Return_pipe_m_flow']] * (row[v_names_loss['Inflow_pipe_t_in']] - row[v_names_loss['Return_pipe_t_out']]) * 4.184, axis = 1) # in kWh
            losses['massflow'] = losses[v_names_loss['Return_pipe_m_flow']]

            losses.reset_index(inplace=True)
            losses['weights'] = losses['s'].map(scenarios_dict.get).apply(lambda x: x)
            total_losses =  losses.filter(regex='tot_losses|weights')[losses['tot_losses'] < 0].prod(axis = 1).sum()
            total_gains =   losses.filter(regex='tot_losses|weights')[losses['tot_losses'] > 0].prod(axis = 1).sum()
            total_generated =  losses.filter(regex='tot_generated|weights').prod(axis = 1).sum()
            total_generated_pos =  losses.filter(regex='tot_generated_source|weights')[losses['tot_generated_source'] > 0].prod(axis = 1).sum()
            total_generated_neg =  losses.filter(regex='tot_generated_source|weights')[losses['tot_generated_source'] < 0].prod(axis = 1).sum()                      
            max_generated = losses['tot_generated'].max() ## maximum heat input required in HX
            max_generated_source = losses['tot_generated_source'].max() ## maximum heat input required in HX
            min_generated_source = losses['tot_generated_source'].min() ## maximum heat input required in HX
            average_mass_flow = losses.filter(regex='Return_pipe_m_flow|weights').prod(axis = 1).sum() / losses['weights'].sum()


            max_generated = losses['tot_generated'].max() ## maximum heat input required in HX
            ## COP WSHP/CENHP/ASHP
            regex_COP = 'COP_h|COP_c'
            COP = operation.filter(regex = regex_COP).copy()
            COP_col = COP.columns.tolist()
            variables_COP = ['COP_h','COP_c']
            v_names_COP = {key:None for key in variables_COP}
            for x in variables_COP:
                index = [i for i, s in enumerate(COP_col) if x in s]
                if len(index) == 1:
                    v_names_COP[x] = COP_col[index[0]]
                elif len(index) == 0:
                    IndexError('Variables of interest in COP have 0 equal names')
                else:
                    IndexError('Variables of interest in COP have repeated names (2 or more')


            COP.reset_index(inplace = True)
            COP['weights'] = COP['s'].map(scenarios_dict.get).apply(lambda x: x)
            COP_h = COP.filter(regex = 'COP_h|weights')[(COP[v_names_COP['COP_h']] != 0) & (COP['weights'] != 0)].prod(axis=1).sum() / COP['weights'][COP[v_names_COP['COP_h']] != 0].sum()
            COP_c = COP.filter(regex = 'COP_c|weights')[(COP[v_names_COP['COP_c']] != 0) & (COP['weights'] != 0)].prod(axis=1).sum() / COP['weights'][COP[v_names_COP['COP_c']] != 0].sum()

        #generation
        #Heat Exchanger addition
        #transmission
        # inflow and return pipe
        op_costs_comp.Pipes = op_costs.filter(regex='Source_el_ind_price|Inflow_pipe_P_pump|weights').prod(axis=1) + op_costs.filter(regex='Source_el_ind_price|Inflow_pipe_P_pump|weights').prod(axis=1)
        
        #Consumption
        # cooling and/or heating
        op_costs_comp.WSHP_LC = (op_costs.filter(regex='Source_el_cons_price|WSHP_LC_C_cons_P_HP_c|weights').prod(axis=1) + op_costs.filter(regex='Source_el_cons_price|WSHP_LC_C_cons_P_HP_h|weights').prod(axis=1))* cons_n
        # heating only
        op_costs_comp.HS_el = (op_costs.filter(regex='Source_el_cons_price|HS_el_C_cons_q_dot|weights').prod(axis=1)) * cons_n
        op_costs_comp.HS_gas = (op_costs.filter(regex='Source_gas_price|Source_gas_use|weights').prod(axis=1)) * cons_n

        op_costs_comp_total = op_costs_comp.sum(axis = 0, skipna = True)
 
        ## Investment & Maintenance costs
        ### Cost per Component
        cost_inv = params.filter(regex = '_p_main|_p_fix|_p_spec').iloc[0]
        cost_list = pd.DataFrame(index = ['Maintenance', 'Investment', 'Operation'], columns = comp_list)

        for comp in comp_list:
            if comp == 'Pipes':
                cost_list.at['Investment',comp] = params.filter(regex = 'Inflow_pipe_pipe_cost').iloc[0,0]*params.filter(regex = 'Inflow_pipe_length').iloc[0,0]
                cost_list.at['Maintenance',comp] = 0
                cost_list.at['Operation',comp] = op_costs_comp_total[comp]

            elif comp == 'Cen_pump':
                cost_list.at['Investment',comp] = 0 * design.filter(regex = 'Cen_pump_Qdot_design').prod(axis=1).iloc[0] * params.filter(regex = 'Cen_pump_p_spec').iloc[0,0] + params.filter(regex = 'Cen_pump_p_fix').iloc[0,0]
                cost_list.at['Maintenance',comp] = 0 * cost_list.at['Investment',comp]/1.3175*0.025 
                cost_list.at['Operation',comp] = 0 * op_costs_comp_total[comp]
            
            elif comp == 'Cen_HX': # Before adding the component to the model
                cost_list.at['Investment',comp] = 115 * max_generated #considers annuity factor for u=25
                cost_list.at['Maintenance',comp] = 0.51 * max_generated # price per kw * needed load
                cost_list.at['Operation',comp] = 0      

            elif comp == 'WSHP_LC':
                cost_list.at['Investment',comp] = (design.filter(regex = comp).filter(regex='WSHP_LC_C_cons_Qdot_design').prod(axis=1).iloc[0] * cost_inv.filter(regex = comp).filter(regex = '_p_spec').iloc[0] +  cost_inv.filter(regex = comp).filter(regex = '_p_fix').iloc[0]) *design.filter(regex = comp).filter(regex='_b_build').iloc[0][0] *cons_n
                cost_list.at['Maintenance',comp] = cost_list.at['Investment',comp]/1.3175*0.025 ## 2.5% removing annuity factor
                cost_list.at['Operation',comp] = op_costs_comp_total[comp]

            else:
                try:
                    cost_list.at['Investment',comp] = (design.filter(regex = comp).filter(regex='_C_cons_Qdot_max').prod(axis=1).iloc[0] * cost_inv.filter(regex = comp).filter(regex = '_p_spec').iloc[0] +  cost_inv.filter(regex = comp).filter(regex = '_p_fix').iloc[0]) * design.filter(regex = comp).filter(regex='_b_build').iloc[0][0]* cons_n
                    cost_list.at['Maintenance',comp] = design.filter(regex = comp).filter(regex='_b_build').iloc[0][0]  * cost_inv.filter(regex = comp).filter(regex = '_p_main').iloc[0]*cons_n
                    cost_list.at['Operation',comp] = op_costs_comp_total[comp]
                except:
                    pass
        
        
        ### Cost per Area
         
        sub_systems = ['Generation','Transmission','Consumers']

        sub_systems_dic = {}

        ### Generation Costs   
        comp_list_gen = ['Cen_HX']
        ### Transmission Costs
        comp_list_trans = ['Pipes']
        ### Consumer Costs
        comp_list_cons = ['HS_gas','HS_el','WSHP_LC'] ## update LC to HX eventually 

        for sub in sub_systems:
            if sub == 'Generation':
                sub_systems_dic['Generation_Investment'] = annuity(cost_list.loc['Investment', comp_list_gen].sum(axis=0),n=30,wacc = 0.03) 
                sub_systems_dic['Generation_Maintenance'] = cost_list.loc['Maintenance', comp_list_gen].sum(axis=0)
                sub_systems_dic['Generation_Operation'] = cost_list.loc['Operation', comp_list_gen].sum(axis=0)      
            if sub == 'Transmission':
                sub_systems_dic['Transmission_Investment'] = annuity(cost_list.loc['Investment', comp_list_trans].sum(axis=0),n=30,wacc = 0.03) 
                sub_systems_dic['Transmission_Maintenance'] = cost_list.loc['Maintenance', comp_list_trans].sum(axis=0)
                sub_systems_dic['Transmission_Operation'] = cost_list.loc['Operation', comp_list_trans].sum(axis=0)       
            if sub == 'Consumers':
                sub_systems_dic['Consumers_Investment'] = annuity(cost_list.loc['Investment', comp_list_cons].sum(axis=0),n=30,wacc = 0.03) 
                sub_systems_dic['Consumers_Maintenance'] = cost_list.loc['Maintenance', comp_list_cons].sum(axis=0)
                sub_systems_dic['Consumers_Operation'] = cost_list.loc['Operation', comp_list_cons].sum(axis=0)

         ### Total costs

        total_costs = {'Investment': annuity(cost_list.loc['Investment'].sum(axis=0),n=30,wacc = 0.03),'Maintenance': cost_list.loc['Maintenance'].sum(axis=0),'Operation':total_commodities}
        Dic_costs[pset_id] = {'Investment_costs': total_costs['Investment'], 'Maintenance_costs': total_costs['Maintenance'], 'Operation_costs': total_costs['Operation'],
        #updated TAC/LCOE
        'upd_tac': sum(total_costs.values()),
        'upd_LCOE': sum(total_costs.values())/tot_demand,
        #Demands
        'tot_dem': tot_demand,
        'h_dem':heat_tot_demand,
        'c_dem':cool_tot_demand,
        'max_generated':max_generated_source,
        'min_generated':min_generated_source,
        'tot_generated_source_pos': total_generated_pos,
        'total_generated_source_neg': total_generated_neg,

        #Losses
        'total_losses': total_losses,
        'total_gains' :  total_gains,
        'total_generated': total_generated,

        #Massflow
        'average_pipe_massflow': average_mass_flow,
        #COP
        'COP_h' : COP_h,
        'COP_c' : COP_c,
        'comp_costs': cost_list.to_dict()} | sub_systems_dic
    
    df_costs = pd.DataFrame.from_dict(Dic_costs).T.reset_index(drop = True)


    pset_post = pd.concat([
                        pset.reset_index(drop=True), # original pset
                        df_costs.drop(['comp_costs'],axis = 1).reset_index(drop = True), ## Total Costs
                        ], axis = 1)


    col_components = ['Generation_Investment','Generation_Maintenance','Generation_Operation',
                'Transmission_Investment','Transmission_Maintenance','Transmission_Operation',
                'Consumers_Investment','Consumers_Maintenance','Consumers_Operation',
                ]

    pickledir_post = parent_path / '_results' / str(run + '_post_processed' + '.pk')

    with pickledir_post.open('wb') as fp:
        pickle.dump((pset_post), fp)

    return pset_post

def multi_solve(pset,run_label,timelim_inp):
    import gurobipy as py
    from pyomo.environ import value
    from comando.interfaces.gurobi import to_gurobi
    from _scripts.model import instantiate_model, run_model
    
    #### Solving options
    #### Timelimit

    if timelim_inp != None:
        timelim_val = timelim_inp
    else:
        timelim_val = 1200
    options = dict(  # Options assuming Gurobi 10.0.1
                Seed=SEED,
                NonConvex=2,
                MIPGap=0.01,
                MIPFocus=2,
                timelimit=timelim_val,
                NumericFocus=0,
                OutputFlag=0,
    )

    ES,names = instantiate_model(label = run_label, ashp_mode = 'cool_only', dt_max = pset['dt_max'], n_seg=10) ## if n_seg is modified, it also has to been added to the input of the model
    if pset['typical_days'] == None: # to pass None argument to run_model function)
        typical_days = None
    else:
        typical_days = pset['typical_days']
    ### Adding ES object to pset_dict
    P = run_model(ES,operational_calc=False, case = pset['region'], typ_days = typical_days, interface = 'gurobi', run_label = run_label, **pset) ## Remember to instantiate model
    m = to_gurobi(P)
    env = py.Env()
    m.env = env
    demands = {}
    details = {}
    try:
        m.solve(**options)
        tac = value(m.objVal)
        status = m.Status
        h_demands = sum(p for p in P.parameters if p.name.endswith('heating_Qdot'))
        c_demands = sum(p for p in P.parameters if p.name.endswith('cooling_Qdot'))
        demands['heating'] = P.weighted_sum(h_demands, symbolic=False) / 1e3
        demands['cooling'] = P.weighted_sum(c_demands, symbolic=False) / 1e3
        demands['total'] = P.weighted_sum(h_demands+c_demands, symbolic=False) / 1e3
        details['design'] = P.design.transpose()
        details['operation'] = P.operation.transpose()
        details['parameters'] =  P.data.getter()
        try:
            details['scenarios'] = P.scenario_weights.to_frame()
        except:
            pass
    except:
        tac = 0
        status = 3
    
    env = None
    m = None
    ES = None
    P = None
    return tac, status, demands, details

def multi_processing(pset,run,**kwargs):
    import pickle
    import pandas as pd
    import numpy as np
    from multiprocessing import Pool,cpu_count
    import gurobipy as py
    from itertools import repeat

    ### Check if only runs that reached timelimit are rerun
    timelim = kwargs.get('timelim', None)

    if timelim == None:
        pickledir_post = parent_path / '_results' / str(run + '_processed.pk')
        timelim_inp = timelim
    elif type(timelim) is not str: ## In case one or more model(s) didnt reach optimality in the first set of runs.
        ori_pset = pset.copy(deep = True)
        timelim_inp = timelim
        pset = pset[(pset['status_list'] == 'Time Limit')] ### could be used for runs that were defined as unfeasible by the solver, requiring a run with a higher NumericFocus parameter or different conditions.
        print(pset)
        pickledir_post = parent_path / '_results' / str(run + '_REprocessed.pk')
    else:
        ValueError('timelim must be None or a number')

    ES_label_list = ['ES_' + str(x) + '_' for x in pset.index]
    pset_dict_list = [x for x in pset['pset_dict']]
    pset_id_list = [x for x in pset['_pset_id']]
    new_calc_dirs_list = [str(parent_path) + x.split(local_folder)[1] for x in pset['_calc_dir']]
    #### empty list for variables
    LCOE_list = []
    str_status_list = []
    unfeasible_list = []
    h_demand_list = []
    c_demand_list = []
    tot_demand_list = []

    #### Optimisation run
    n_cores = cpu_count()
    with Pool(processes = n_cores,maxtasksperchild=10) as pool:
        tac_list, status_list, demand_list, details_list = zip(*pool.starmap(multi_solve, zip(pset_dict_list,ES_label_list,repeat(timelim_inp)),chunksize = 100))

    #### from tuple to list
    tac_list = [*tac_list]
    status_list = [*status_list]
    demand_list = [*demand_list]
    details_list = [*details_list]

    #### Data processing and storage
    for c_dir,pset_id,pset_it,demand,details,tac,status in zip(new_calc_dirs_list,pset_id_list,pset_dict_list,demand_list,details_list,tac_list,status_list):

        targetdir = c_dir + '/' + pset_id
        length_string = 'length = %4.1f m' %pset_it['length']
        ncon_string = 'ncon = %5d' %pset_it['n_cons']
        Tin_string ='Tin = %5.2f K' %pset_it['T_source']
        DN_string = 'DN = %5.2f mm' %pset_it['DN']    

        if status == 2: #optimal solution found - for given MIPgap
            str_status = 'Optimal'
            tot_demand = demand['total']
            h_demand_list.append(demand['heating'])
            c_demand_list.append(demand['cooling'])
            tot_demand_list.append(tot_demand)
            details['design'].to_csv(targetdir + '/' + 'design_var.csv')
            details['operation'].to_csv(targetdir + '/' + 'operation_var.csv')
            details['parameters'].to_csv(targetdir + '/' + 'parameters_inp.csv')

        elif status == 3:
            str_status = 'Unfeasible'
            unfeasible_list.append(f'Unfeasible solution for case: {length_string} | {ncon_string} | {Tin_string} | {DN_string}')
            tot_demand = 1
            h_demand_list.append(1)
            c_demand_list.append(1)
            tot_demand_list.append(tot_demand)

        else:
            str_status = 'Time Limit'
            unfeasible_list.append(f'Timelimit for case: {length_string} | {ncon_string} | {Tin_string} | {DN_string}, status = {status}, pset id = {pset_id}')
            tot_demand = demand['total']
            h_demand_list.append(demand['heating'])
            c_demand_list.append(demand['cooling'])
            tot_demand_list.append(tot_demand)
            details['design'].to_csv(targetdir + '/' + 'design_var.csv')
            details['operation'].to_csv(targetdir + '/' + 'operation_var.csv')
            details['parameters'].to_csv(targetdir + '/' + 'parameters_inp.csv')
 
        LCOE_list.append(tac/tot_demand)
        str_status_list.append(str_status)

        with (Path(targetdir + '/' + 'details.pk')).open('wb') as fp:
            pickle.dump((details), fp)


    pset['tac'] = tac_list
    pset['LCOE'] = LCOE_list
    pset['status_list'] = str_status_list
    pset['tot_dem'] = tot_demand_list
    pset['h_dem'] = h_demand_list
    pset['c_dem'] = c_demand_list

    if timelim == None:
        with open(c_dir + '/' + 'Unfeasible.txt', 'w') as f:
            for line in unfeasible_list:
                    f.write(f"{line}\n")
        out_pset = pset

    else: ## timelim not being a string has already been checked
        with open(c_dir + '/' + 'checklog_timelimit.txt', 'w') as f:
            for line in unfeasible_list:
                    f.write(f"{line}\n")
        ori_pset.set_index('_pset_id', drop = False)
        ori_pset.update(pset)
        ori_pset.reset_index(drop = True, inplace = True)
        out_pset = ori_pset

    with pickledir_post.open('wb') as fp:
        pickle.dump((out_pset), fp)
    
    return out_pset

########
### Model Builder & execution
########    

if __name__ == '__main__' and reprocessing == False:
    import psweep as ps
    from _scripts.model import instantiate_model, run_model
    from pathlib import Path
    import pandas as pd
    import os

    ######################
    # Commodities and Piping Sizing
    ######################    

    elec_ind_NRW = 18.33
    elec_cons_NRW = 32.79
    elec_fac_NRW = elec_ind_NRW/elec_cons_NRW
    gas_cons_NRW = 8.06
    gas_ind_NRW = 5.03
    gas_fac_NRW = gas_ind_NRW/gas_cons_NRW

    elec_ind_VIC = 27.09/1.5
    elec_cons_VIC = 33.55/1.5
    elec_fac_VIC = elec_ind_VIC/elec_cons_VIC
    gas_cons_VIC = 10.24/1.5
    gas_ind_VIC = 7.70/1.5
    gas_fac_VIC = gas_ind_VIC/gas_cons_VIC

    

    data_path = parent_path / '_data'
    data_path_VIC = data_path / 'VIC'
    data_path_NRW = data_path / 'NRW'

    #Pipe dimensioning

    if ins == 'NoIns':
        pipe_data = pd.read_csv(data_path / 'pipe_PE_prop.csv', sep=',') # No insulation PE pipes
    elif ins == 'S1' or ins == 'S2' or ins == 'S3':
        pipe_data = pd.read_csv(data_path / 'pipe_steel_prop.csv', sep=',') # Steel pipes with different insulation thickness

    pipe_data_VIC = pipe_data
    q_c_vic = 6.43 ## maximum cooling load
    q_h_vic = 7.29 ## maximum heating load

    pipe_data_NRW = pipe_data
    q_c_nrw = 1.68 ## maximum cooling load
    q_h_nrw = 5.97 ## maximum heating load

    # Region Checker
    if region == 'VIC':
        q_c = q_c_vic
        q_h = q_h_vic
        pipe_data = pipe_data_VIC
        elec_cons = elec_cons_VIC
        elec_ind = elec_ind_VIC
        elec_fac = elec_fac_VIC
        gas_cons = gas_cons_VIC
        gas_ind = gas_ind_VIC
        gas_fac = gas_fac_VIC
        elec_gas_ratio = elec_cons / gas_cons

    elif region =='NRW':
        q_c = q_c_nrw
        q_h = q_h_nrw
        pipe_data = pipe_data_NRW
        elec_cons = elec_cons_NRW
        elec_ind = elec_ind_NRW
        elec_fac = elec_fac_NRW
        gas_cons = gas_cons_NRW
        gas_ind = gas_ind_NRW
        gas_fac = gas_fac_NRW
        elec_gas_ratio = elec_cons / gas_cons

    ######################
    # Energy System Parameters
    ######################

    #list generation
    DN_long =  []
    TmaxoutC_long = []
    mmax_long = []
    dt_max_long = []
    T_long = []
    n_long = []
    L_long = []
    label_vect = [str(l) for l in len_vect]
    for Tin in Tin_vect:
        for n in n_con:
            for L in len_vect:
                T_long.append(Tin)
                n_long.append(n)
                L_long.append(L)
                dt_max = 25
                dt_max_long.append(dt_max)
                DN_it,TmaxoutC,mmax = DN_calc(dt_max,Tin,q_c,q_h,n,pipe_data,L,'yes')
                TmaxoutC_long.append(TmaxoutC)
                DN_long.append(DN_it)
                mmax_long.append(mmax)

    if append_base == True:
        el_price_list.append(elec_cons/100)
    gas_price_list = [gas_cons/100]
    ins_param = ps.plist('ins',[ins])
    dt_max_param = ps.plist('dt_max',dt_max_long)
    el_price_param = ps.plist('el_price',el_price_list)
    el_fac_param = ps.plist('elec_factor', [elec_fac])
    gas_price_param = ps.plist('gas_price',gas_price_list)
    len_param = ps.plist('length',L_long)
    n_con_param = ps.plist('n_cons',n_long)
    Tin_param =  ps.plist('T_source',T_long)
    DN_temp_param=ps.plist('DN',DN_long)


    ######################
    # Generate Model Instance & Run
    ######################


    typical_days = ps.plist('typical_days',[typ_day_file])  # Either 'None' for a complete timeseries (8760 hrs) or the filename with the aggregation
    region_param = ps.plist('region', [region])
        ### running multiprocessing version
    params = ps.pgrid(el_price_param,gas_price_param,zip(len_param,n_con_param,DN_temp_param,Tin_param,dt_max_param),el_fac_param,typical_days,region_param,ins_param)
    df = ps.run_local(calc_func_mult, 
                    params = params,
                    database_dir = str(parent_path / '_results'),
                    calc_dir= str(parent_path / '_results' / run),
                    database_basename= run + '.pk',
                    )
    print(df)
    df = multi_processing(df,run)
    print(df)
    df_post = post_process_cost(df,run)
    print(df_post)


########
### Reprocessing -- if reprocessing / repostprocessing of the models is required
########   

if __name__ == '__main__' and reprocessing == True and timelimit == False:

    import psweep as ps
    from pathlib import Path
    import pandas as pd
    import os  
    import pickle
    from _scripts.model import instantiate_model, run_model

    print('Multiprocessing ' + str(' '.join(repro_name_list)) + ' Results. Cancel if it is not the intended script execution.')
    #reprocessing already ran sets
    results_path = parent_path / '_results'

    os.chdir(str(parent_path))

    ##### list of files
    if not repro_name_list:
        pass
    else:
        for name in repro_name_list:
            dirfile = [f for f in results_path.glob('*' + name + '_processed' +'.pk')]
            with open(f'{dirfile[0]}', 'rb') as f: # last generated pickle
                db = pickle.load(f)
                db = multi_processing(db,name)
                print(db)
                db_post = post_process_cost(db,name)
                print(db)

########
### Reprocessing - Timelimit
########

if __name__ == '__main__' and reprocessing == True and timelimit == True:

    import psweep as ps
    from pathlib import Path
    import pandas as pd
    import os  
    import pickle
    from _scripts.model import instantiate_model, run_model

    timelimit = 600
    print('Multiprocessing ' + str(' '.join(repro_timelimit_list)) + ' Results for runs that reached timelimit, setting value to timelimit = ' + str(timelimit) + 's')
    #reprocessing already ran sets
    results_path = parent_path / '_results'

    os.chdir(str(parent_path))

    ##### list of files
    if not repro_timelimit_list:
        pass
    else:
        for name in repro_timelimit_list:
            dirfile = [f for f in results_path.glob('*' + name +'_processed.pk')]
            with open(f'{dirfile[0]}', 'rb') as f: # last generated pickle
                db = pickle.load(f)
                db = multi_processing(db,name,timelim = timelimit)
                print(db)
                db_post = post_process_cost(db,name)
                print(db)