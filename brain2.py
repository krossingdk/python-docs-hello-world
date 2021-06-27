# Copyright 2010-2018 Google LLC
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# Copyright 2010-2018 Google LLC
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
"""Creates a shift scheduling problem and solves it."""

from __future__ import print_function
from itertools import combinations
import argparse
import math

from google.protobuf import text_format
from ortools.sat.python import cp_model
import sys
sys.stdout = open('brainsolution.txt','wt')

import json


def test_overlap(t1_st, t1_end, t2_st, t2_end):

    def convert_to_minutes(t_str):
        hours, minutes = t_str.split(':')
        return 60*int(hours)+int(minutes)

    t1_st = convert_to_minutes(t1_st)
    t1_end = convert_to_minutes(t1_end)
    t2_st = convert_to_minutes(t2_st)
    t2_end = convert_to_minutes(t2_end)

    # Check for wrapping time differences
    if t1_end < t1_st:
        if t2_end < t2_st:
        # Both wrap, therefore they overlap at midnight
            return True
        # t2 doesn't wrap. Therefore t1 has to start after t2 and end before
        return t1_st < t2_end or t2_st < t1_end

    if t2_end < t2_st:
        # only t2 wraps. Same as before, just reversed
        return t2_st < t1_end or t1_st < t2_end

    # They don't wrap and the start of one comes after the end of the other,
    # therefore they don't overlap
    if t1_st >= t2_end or t2_st >= t1_end:
        return False
    # In all other cases, they have to overlap
    return True

def in_nested_list(my_list, item):
    """
    Determines if an item is in my_list, even if nested in a lower-level list.
    """
    if item in my_list:
        return True
    else:
        return any(in_nested_list(sublist, item) for sublist in my_list if isinstance(sublist, list))


PARSER = argparse.ArgumentParser()
PARSER.add_argument(
    '--output_proto',
    default="",
    help='Output file to write the cp_model'
    'proto to.')
PARSER.add_argument('--params', default="", help='Sat solver parameters.')


def negated_bounded_span(works, start, length):
    """Filters an isolated sub-sequence of variables assined to True.

  Extract the span of Boolean variables [start, start + length), negate them,
  and if there is variables to the left/right of this span, surround the span by
  them in non negated form.

  Args:
    works: a list of variables to extract the span from.
    start: the start to the span.
    length: the length of the span.

  Returns:
    a list of variables which conjunction will be false if the sub-list is
    assigned to True, and correctly bounded by variables assigned to False,
    or by the start or end of works.
  """
    sequence = []
    # Left border (start of works, or works[start - 1])
    if start > 0:
        sequence.append(works[start - 1])
    for i in range(length):
        sequence.append(works[start + i].Not())
    # Right border (end of works or works[start + length])
    if start + length < len(works):
        sequence.append(works[start + length])
    return sequence


def add_soft_sequence_constraint(model, works, hard_min, soft_min, min_cost,
                                 soft_max, hard_max, max_cost, prefix):
    """Sequence constraint on true variables with soft and hard bounds.

  This constraint look at every maximal contiguous sequence of variables
  assigned to true. If forbids sequence of length < hard_min or > hard_max.
  Then it creates penalty terms if the length is < soft_min or > soft_max.

  Args:
    model: the sequence constraint is built on this model.
    works: a list of Boolean variables.
    hard_min: any sequence of true variables must have a length of at least
      hard_min.
    soft_min: any sequence should have a length of at least soft_min, or a
      linear penalty on the delta will be added to the objective.
    min_cost: the coefficient of the linear penalty if the length is less than
      soft_min.
    soft_max: any sequence should have a length of at most soft_max, or a linear
      penalty on the delta will be added to the objective.
    hard_max: any sequence of true variables must have a length of at most
      hard_max.
    max_cost: the coefficient of the linear penalty if the length is more than
      soft_max.
    prefix: a base name for penalty literals.

  Returns:
    a tuple (variables_list, coefficient_list) containing the different
    penalties created by the sequence constraint.
  """
    cost_literals = []
    cost_coefficients = []

    # Forbid sequences that are too short.
    for length in range(1, hard_min):
        for start in range(len(works) - length + 1):
            model.AddBoolOr(negated_bounded_span(works, start, length))

    # Penalize sequences that are below the soft limit.
    if min_cost > 0:
        for length in range(hard_min, soft_min):
            for start in range(len(works) - length + 1):
                span = negated_bounded_span(works, start, length)
                name = ': under_span(start=%i, length=%i)' % (start, length)
                lit = model.NewBoolVar(prefix + name)
                span.append(lit)
                model.AddBoolOr(span)
                cost_literals.append(lit)
                # We filter exactly the sequence with a short length.
                # The penalty is proportional to the delta with soft_min.
                cost_coefficients.append(min_cost * (soft_min - length))

    # Penalize sequences that are above the soft limit.
    if max_cost > 0:
        for length in range(soft_max + 1, hard_max + 1):
            for start in range(len(works) - length + 1):
                span = negated_bounded_span(works, start, length)
                name = ': over_span(start=%i, length=%i)' % (start, length)
                lit = model.NewBoolVar(prefix + name)
                span.append(lit)
                model.AddBoolOr(span)
                cost_literals.append(lit)
                # Cost paid is max_cost * excess length.
                cost_coefficients.append(max_cost * (length - soft_max))

    # Just forbid any sequence of true variables with length hard_max + 1
    for start in range(len(works) - hard_max):
        model.AddBoolOr(
            [works[i].Not() for i in range(start, start + hard_max + 1)])
    return cost_literals, cost_coefficients


def add_soft_sum_constraint(model, works, hard_min, soft_min, min_cost,
                            soft_max, hard_max, max_cost, prefix):
    """Sum constraint with soft and hard bounds.

  This constraint counts the variables assigned to true from works.
  If forbids sum < hard_min or > hard_max.
  Then it creates penalty terms if the sum is < soft_min or > soft_max.

  Args:
    model: the sequence constraint is built on this model.
    works: a list of Boolean variables.
    hard_min: any sequence of true variables must have a sum of at least
      hard_min.
    soft_min: any sequence should have a sum of at least soft_min, or a linear
      penalty on the delta will be added to the objective.
    min_cost: the coefficient of the linear penalty if the sum is less than
      soft_min.
    soft_max: any sequence should have a sum of at most soft_max, or a linear
      penalty on the delta will be added to the objective.
    hard_max: any sequence of true variables must have a sum of at most
      hard_max.
    max_cost: the coefficient of the linear penalty if the sum is more than
      soft_max.
    prefix: a base name for penalty variables.

  Returns:
    a tuple (variables_list, coefficient_list) containing the different
    penalties created by the sequence constraint.
  """
    cost_variables = []
    cost_coefficients = []
    sum_var = model.NewIntVar(hard_min, hard_max, '')
    # This adds the hard constraints on the sum.
    model.Add(sum_var == sum(works))
    # Penalize sums below the soft_min target.
        
    if soft_min > hard_min and min_cost > 0:
        delta = model.NewIntVar(-len(works), len(works), '')
        model.Add(delta == soft_min - sum_var)
        # TODO(user): Compare efficiency with only excess >= soft_min - sum_var.
        excess = model.NewIntVar(0, 7, prefix + ': under_sum')
        model.AddMaxEquality(excess, [delta, 0])
        cost_variables.append(excess)
        cost_coefficients.append(min_cost)

    # Penalize sums above the soft_max target.
    if soft_max < hard_max and max_cost > 0:
        delta = model.NewIntVar(-7, 7, '')
        model.Add(delta == sum_var - soft_max)
        excess = model.NewIntVar(0, 7, prefix + ': over_sum')
        model.AddMaxEquality(excess, [delta, 0])
        cost_variables.append(excess)
        cost_coefficients.append(max_cost)

    return cost_variables, cost_coefficients


def solve_shift_scheduling(params, output_proto):
    import json
    import pickle
   
    with open('data.json', 'r') as f:
        data=json.load(f)
   
    skills = json.loads(data['skills'])
    requests = json.loads(data['requests'])
    cover_demands = json.loads(data['cover_demands'])
    weekly_sum_constraints = json.loads(data['weekly_sum_constraints'])
    funktioner_id_all = json.loads(data['funktioner_id_all'])
    print('funktioner_id_all',funktioner_id_all)
    funktioner_funktion_all = json.loads(data['funktioner_funktion_all'])
    funktioner_liste = list(range(0,len(funktioner_funktion_all)))
    print('funktioner_liste', funktioner_liste)
    funktioner_dag_index_sleep = json.loads(data['funktioner_dag_index_sleep'])
    print('valgte dagfunktioner', funktioner_dag_index_sleep)
    funktioner_vagt_index_sleep = json.loads(data['funktioner_vagt_index_sleep'])
    print('valgte vagter', funktioner_vagt_index_sleep )
    funktioner_vagt_dag_index_sleep = json.loads(data['funktioner_vagt_dag_index_sleep'])
    print('valgte dagvagter', funktioner_vagt_index_sleep )
    funktioner_vagt_night_index_sleep = json.loads(data['funktioner_vagt_night_index_sleep'])
    print('valgte nattevagter',funktioner_vagt_night_index_sleep)
    night_shifts = funktioner_vagt_night_index_sleep
    #daytime_oncall_shifts = funktioner_vagt_dag_index_sleep
    print ('nattevagter' , night_shifts)
    print('antal nattevagter', len(night_shifts))
    
    equalfunctions = json.loads(data['equalfunctions'])
    equalfunctionsansatte = json.loads(data['equalfunctionsansatte'])
    
    print('equalfunntionansatte', equalfunctionsansatte)
    #timer (ansat:arbejdstimer(fri ønkser+timer for ikke valgte funktioner), dagvagter, nattevagter, resttimer
    timer = json.loads(data['timer'])
    print('timer -arb timer dagvagter, nattevagter, rest', timer)
    calc_time = int(data['calc_time'])
    #calc_time=300
    print('tænketid: ',calc_time)
    #print('calc_time', data['calc_time'])
    #print ('data', data)
    
    ansatte = json.loads(data['ansatte'])
    offset = json.loads(data['offset'])
    period_length = data['period_length']
    number_of_weeks = data['number_of_weeks']
    fixed_assignments = json.loads(data['fixed_assignments'])
    desired_shift_transitions = json.loads(data['desired_shift_transitions'])
    desired_day_transistions = json.loads(data['desired_day_transistions'])
    ansat_arbejdstid = json.loads(data['ansat_arbejdstid'])
    ugedag_nummer = json.loads(data['ugedag_nummer'])
    ugedag_dag = json.loads(data['ugedag_dag'])
    dates = json.loads(data['dates'])
    funktioner_tider = json.loads(data['funktioner_tider'])
    funktioner_varighed = json.loads(data['funktioner_varighed'])
    normtid= json.loads(data['normtid'])

   
    funktioner_max_uge = json.loads(data['funktioner_max_uge'])
    funktioner_allow_overlap = json.loads(data['allow_overlap'])
    print('funktioner_allow_overlap', funktioner_allow_overlap);
    
    funktioner_no_overlap= [shift for shift in funktioner_liste if shift not in funktioner_allow_overlap]
    print('funktioner_no_overlap',funktioner_no_overlap)
    
    funktioner_allow_combination = json.loads(data['allow_combination'])
    funktioner_no_combination= [shift for shift in funktioner_liste if shift not in funktioner_allow_combination]
    print('funktioner_no_combination',funktioner_no_combination)
   
    print('normtid', normtid)
   
   
    
    fixed_assignments_shift_day = []
    for x in fixed_assignments:
        fixed_assignments_shift_day.append([x[1],x[2]])
    
  
    print ('fixed', fixed_assignments)    
    print('funktioner liste', funktioner_liste)
    print('equalfunctions', equalfunctions)
    print('funktionstider', funktioner_tider)
    print('funktionsvarighed', funktioner_varighed)
    print('funktioner max', funktioner_max_uge)
    print('dv var', funktioner_varighed[str(1)])
   
    print('fkt tid 0', funktioner_tider['1'][2])
    print('dates', dates, 'funktioner tider' , funktioner_tider, 'ugedag_nummer',  ugedag_nummer)
    print('desired_day_transistions', desired_day_transistions)
    print(desired_shift_transitions, 'desired_shift_transitions')
    print ('offset', offset)
    print ('fixed assignments', fixed_assignments)
    print('skills', skills)
    print ('funktioner id', funktioner_id_all)
    print ('arbejdstid', '1', int(math.floor(ansat_arbejdstid[0][1])*100/37))
    print ('ansatte', ansatte)
    
    
    print('antal ansatte', len(ansatte))
    print(weekly_sum_constraints)
    print('alle funktioner', funktioner_funktion_all)
    print ('cover_demands', cover_demands)
    print('requests', requests)
    
    
    
   
        
    
   
    # Data
        # Shift constraints on continuous sequence :
    #     (shift, hard_min, soft_min, min_penalty,
    #             soft_max, hard_max, max_penalty)
    shift_constraints = [
        # One or two consecutive days of rest, this is a hard constraint.
        #(0, 1, 2, 2, 5, 7, 0),
        # betweem 2 and 3 consecutive days of night shifts, 1 and 4 are
        # possible but penalized.
     #   (2, 0, 1, 2, 0,1, 0),
     #  (1, 0, 1, 2, 0,1, 0),
     #   (3, 0, 1, 2, 0,1, 0),
      #  (7, 0, 1, 2, 0,1, 0),
    ]
    
    #first_shift,weekday_first, second_shift, weekday _second
    #desired_transitions = []
     #  (2, 4, 1, 6 ),
    #     (1, 5, 2, 6 ),
    #    (5, 4, 5, 5 ),
    #    (5, 5, 5, 6 ),
    #    
    #    (7, 4, 3, 6 ),
    #    (3, 5, 7, 6 ),
        
         
   # ]
#wanted    
    rewarded_transitions = []
        # B_DV skal efterfølges af B_NV straffes 4
 #       (3, 7, 4),
        # F_DV skal efterfølges af F_NV 0=hard contraint
     #   (1, 2, 4),
      
        # Night to morning is forbidden.
        #(3, 1, 0),
   # ]
   
    
    ansatte_dummy = ansatte.copy()
    ansatte_dummy.append('?')
    print('ansatte', ansatte, 'dummy', ansatte_dummy, 'længde', len(ansatte_dummy))
    employees= ansatte.copy()
    employees_dummy=ansatte_dummy.copy()
    num_employees = len(employees)
    num_employees_dummy = len( employees_dummy)
    print('ansatte dummy', ansatte_dummy)
    print('empl dummy ', employees_dummy)
    num_weeks = number_of_weeks
    print('num_weeks',num_weeks)
    shifts = funktioner_funktion_all
    num_shifts = len(shifts)
    print('antal funktioner', num_shifts, 'alle valte funktioner', funktioner_funktion_all)
    
    #shifts = funktioner_funktion_all
    #num_shifts = len(shifts)
    #print('antal funktioner', num_shifts)
    
    num_days = int(period_length)
    print('antal dage', int(period_length) )
   
    
    
# """Solves the shift scheduling problem.""" 
    model = cp_model.CpModel()

    work = {}
    for e in range(num_employees_dummy):
        for s in range(num_shifts):
            for d in range(num_days):
                work[e, s, d] = model.NewBoolVar('work%i_%i_%i' % (e, s, d))
       

    

    # Linear terms of the objective in a minimization context.
    obj_int_vars = []
    obj_int_coeffs = []
    obj_bool_vars = []
    obj_bool_coeffs = []


# Max shift per e per day.
    for e in range(num_employees):
        for d in range(num_days):
            model.Add(sum(work[e, s, d] for s in range(num_shifts)) <= 2)
            model.Add(sum(work[e, s, d] for s in range(num_shifts)) >= 1)
            model.Add(sum(work[e, s, d] for s in funktioner_no_combination) <= 1)
            #model.Add(sum(work[e, s, d] for s in range(funktioner_allow_combination)) <= 2)
            
            #model.Add((sum(work[e, s, d] for s in night_shifts)+work[e,0,d] ) <= 1)
            #if shift duration < 6 hours shift must be combined with another shift
            for s in range(1, num_shifts):
                if (float(funktioner_tider[str(s)][2]) < 6):
                    model.Add(sum(work[e, s, d] for s in range(num_shifts)) > 1).OnlyEnforceIf(work[e, s, d])  
            #model.Add(sum(work[e, s, d]* int(funktioner_varighed[str(s)]) for s in range(num_shifts)) > 6).OnlyEnforceIf(work[e, s, d])  

        
# constraint if one shift duration < 6 hours other shift same day is not free (0)
    for s in range(1, num_shifts):
        if (float(funktioner_tider[str(s)][2]) < 6):
            #print('tider max', funktioner_tider[str(s)])
            for d in range(num_days):
                for e in range(num_employees):
                    model.AddBoolOr([work[e,s,d].Not(), work[e,0,d].Not()])
                                

#dummy
    for d in range(num_days):        
        model.Add(sum(work[num_employees_dummy-1, s, d] for s in range(num_shifts)) >= 1)   

#dummy penalty    
    num_dummy= sum(work[num_employees_dummy-1, s, d] for s in range(1,num_shifts) for d in range(offset, num_days))
    num_dummy_v = model.NewIntVar(0,num_days*num_shifts,'num_dummies')
    model.Add(num_dummy_v==num_dummy)
    obj_int_vars.append(num_dummy_v)
    obj_int_coeffs.append(300)

# constraint where 2 shifts are assigned to employee, these must not overlap
    
    for (s1, s2) in combinations(funktioner_no_overlap, 2):
        if test_overlap(funktioner_tider[str(s1)][0], funktioner_tider[str(s1)][1], funktioner_tider[str(s2)][0], funktioner_tider[str(s2)][1]):
            for d in range(num_days):
                for e in range(num_employees):
                    model.AddBoolOr([work[e,s1,d].Not(), work[e,s2,d].Not()])

# constraint if overlap make sure 2. shift same day is not 0 to avoid e getting overlap when free
    for s in funktioner_allow_overlap:
        for d in range(num_days):
            for e in range(num_employees):
                model.AddBoolOr([work[e,s,d].Not(), work[e,0,d].Not()])                    
                    
                        
    
# Fixed assignments.
    for e, s, d in fixed_assignments:
        #if d>offset:
        if d>(offset-1):
            model.Add(work[e, s, d] == 1)
    

# Employee requests
    for e, s, d, w in requests:
        if d>(offset-1):
        #if d>offset:    
            if w==0:
                #print(ansatte[e],s,d,w, 'tider', funktioner_tider[str(s)][0], funktioner_tider[str(s)][1])
                model.Add(work[e,s,d] == 1)
                
               
            else:
                obj_bool_vars.append(work[e, s, d])
                obj_bool_coeffs.append(w)

   
#constraint workers needed = workers assigned
    for d in range (num_days):
        for s in range(1, num_shifts):
            #if [s,d] not in fixed_assignments_shift_day: #avoid fixed assignment conflict
                needed = cover_demands[d][s-1]
                model.Add(sum(work[e, s, d] for e in range(num_employees_dummy)) == needed )  

# rest after nightshift
    for night_shift in night_shifts:
        for e in range(num_employees):
            for d in range(num_days - 1):
             #   for s in range(len(shifts)):
              #      if s != night_shift:
                        #sleep_transition = [
                         #   work[e, night_shift, d].Not(),
                          #  work[e, 0, d + 1].Not()
                        #]
                       #model.AddBoolOr(sleep_transition)
                       model.Add(work[e, 0, d + 1]==1).OnlyEnforceIf(work[e, night_shift, d])
                       #model.AddBoolOr([work[e, 0, d].Not(), work[e, night_shift, d].Not()])
                       
   #  model.Add(work[e, 0, d]==1 for d in range(num_days) for e in range(num_empl) ).OnlyEnforceIf(work[e, night_shift, d])                        

#max 1 nightshift indenfor 7 dage
    for e in range(num_employees):
        for d in range(offset, num_days-7):
            model.Add(sum(work[e, s, x] for s in (night_shifts) for x in range(d,d+7)) <= 1)


#    for e in range(num_employees):
#            for s, max in funktioner_max_uge:
#                for d in range(num_days-7):
#                    model.Add(sum(work[e, s, x] for x in range(d,d+7)) <= 7)


#max 1 nightshift indenfor 7 dage
#    for e in range(num_employees):
#        for d in range(num_days-7):
#            nights_worked_sevendays = [work[e, s, x] for s in (night_shifts) for x in range(d,d+7)]
#            max_nights = 1
#            nights_worked = model.NewIntVar(0, 7, '')
#            model.Add(nights_worked == sum(nights_worked_sevendays))
#            over_penalty_nights = 2
#            name = 'excess_nights(empl=%i, shift=%i, day=%i)' % (e, s, d)
#            excess_nights = model.NewIntVar(-1, 7, name)
#            model.Add(excess_nights == nights_worked - max_nights)
#            obj_int_vars.append(excess_nights)
#            obj_int_coeffs.append(over_penalty_nights)

#max shifts
    for e in range(num_employees):
        for s, max in funktioner_max_uge:
            for w in range(num_weeks):
                shifts_worked_sevendays = sum(work[e, s, x] for x in range(w*7,w*7+6) )
                max_shifts = int(max)
                shifts_worked = model.NewIntVar(0, 7, '')
                model.Add(shifts_worked == shifts_worked_sevendays)
                over_penalty_shifts = 4
                name = 'excess_shifts(empl=%i, shift=%i, day=%i)' % (e, s, d)
                excess_shifts = model.NewIntVar(-7, 7, name)
                model.Add(excess_shifts == (max_shifts - shifts_worked))
                obj_int_vars.append(excess_shifts)
                obj_int_coeffs.append(over_penalty_shifts)




#fordel timer, max timer   
    for e in range(num_employees):
        num_hours_worked = sum( work[e, s, d] * int(funktioner_varighed[str(s)]) for s in range(1, num_shifts) for d in range (offset, num_days) )
        max_workhour = int(37*normtid / ansat_arbejdstid[e][1])
        #hours_worked = model.NewIntVar(0,num_days*24,'')
        hours_worked = model.NewIntVar(0, max_workhour,'')
        model.Add(hours_worked==num_hours_worked)
        penalty = -10
        name = 'excess_workhours(empl=%i)' % (e)
        excess_workhours=model.NewIntVar(0,num_days*24, name)
        model.Add(excess_workhours== max_workhour - hours_worked) 
        obj_int_vars.append(excess_workhours)
        obj_int_coeffs.append(penalty)
        
        #norm_hours=int(normtid * ansat_arbejdstid[s][1]/37)
        #penalty = -
        #name = 'less_workhours(empl=%i)' % (e)
        ##less_workhours=model.NewIntVar(0,num_days*24, name)
        ##model.Add(less_workhours==max_workhour - hours_worked) 
        ##obj_int_vars.append(less_workhours)
        ##obj_int_coeffs.append(penalty)
        
            
#penalty for 2 weekends in a row
#    for e in range(num_employees):
#        for w in range(num_weeks-2):
#            num_two_consecutive_weekends = sum(work[e, s, w*7+5]+work[e, s, w*7+6]+work[e, s, w*7+12]+work[e, s, w*7+13] for s in range(1,num_shifts) )
#            num_two_consecutive_weekends_v = model.NewIntVar(0,4,'two_consecutive_weekends_v')
#            model.Add(num_two_consecutive_weekends_v==num_two_consecutive_weekends)
#            obj_int_vars.append(num_two_consecutive_weekends_v)
#            obj_int_coeffs.append(10)
    
    
    
        
#fordel nattevagter via equality for de enkelte mattevagter    
#    sum_of_shifts = {}
#    for e in range(num_employees):
#        for s in night_shifts:
#            sum_of_shifts[(e, s)] = model.NewIntVar(0, num_days, 'sum_of_shifts_%i_%i' % (e, s))
#            model.Add(sum_of_shifts[(e, s)] == sum(work[e, s, d] for d in range(num_days)))
            
#    for s in night_shifts:
#        min_fair_shift = model.NewIntVar(0, num_days, 'min_fair_shift_%i' % s)
#        max_fair_shift = model.NewIntVar(0, num_days, 'max_fair_shift_%i' % s)
#        model.AddMinEquality(min_fair_shift, [sum_of_shifts[(e, s)] for e in range(num_employees)])
#        model.AddMaxEquality(max_fair_shift, [sum_of_shifts[(e, s)] for e in range(num_employees)]) 
#        nightshift_diff = model.NewIntVar(0, num_days, '')
#        model.Add(nightshift_diff==max_fair_shift - min_fair_shift)
#        penalty = 2
#        obj_int_vars.append(nightshift_diff)
#        obj_int_coeffs.append(penalty)


#equalize number of weekends worked inc friday night
    if 2>1:
        ansatte_number_of_weekends_worked_8weeks = json.loads(data['ansatte_number_of_weekends_worked_8weeks'])
        print('ansatte_number_of_weekends_worked_8weeks',ansatte_number_of_weekends_worked_8weeks)
        
        weekends_worked = {}
        weekends = {}
        weekends_in_month = {}
        total_weekends = {}
        trans_avoid_2weekends = {}
        trans_avoid_2weekends_var = {}
        total = {}
        weekenddays_worked_int = {}
        weekenddays_worked = {}
        for e in range(num_employees):
            for w in range(number_of_weeks):
                b=model.NewBoolVar('b')
                weekenddays_worked[e, w]=model.NewBoolVar('weekenddays_worked[%i,%i]' % (e,w))
                weekenddays_worked_int[e, w] = model.NewIntVar(0, number_of_weeks*3,  'weekenddays_worked_int[%i,%i]' % (e,w) )
                weekends[e, w] = model.NewIntVar(0, number_of_weeks,  'weekends[%i,%i]' % (e,w) )
        
                model.Add(weekenddays_worked_int[e, w] ==( sum(work[e, s, w*7+5]+work[e, s, w*7+6] for s in range(1,num_shifts)) +sum(work[e, friday, w*7+4] for friday in night_shifts  )) )
                model.Add(weekenddays_worked_int[e, w]>0).OnlyEnforceIf(b)
                model.Add(weekenddays_worked_int[e, w]==0).OnlyEnforceIf(b.Not())
                
                model.Add(weekenddays_worked[e, w]==1).OnlyEnforceIf(b)
                model.Add(weekenddays_worked[e, w]==0).OnlyEnforceIf(b.Not())
                
                #model.AddDivisionEquality(weekends[e, w],weekenddays_worked_int[e, w],3)
        #avoid two weekends in a row
            if 1<2:
                cost=8
                for wc in range(number_of_weeks-1):
                    trans_avoid_2weekends = [weekenddays_worked[e, wc].Not(), weekenddays_worked[e, wc+1].Not()]
                    trans_avoid_2weekends_var = model.NewBoolVar('trans_avoid_2weekends (employee=%i, week=%i)' % (e, wc))
                    
                    trans_avoid_2weekends.append(trans_avoid_2weekends_var)     
                    model.AddBoolOr(trans_avoid_2weekends)
                    obj_bool_vars.append(trans_avoid_2weekends_var)
                    obj_bool_coeffs.append(cost)
                    
                    
        #max 2 weekends per month
            weekends_in_month[e] = model.NewIntVar(0, 2,  'weekends_in_month[%i]' % e )
            model.Add(weekends_in_month[e] == sum(weekenddays_worked[e, w] for w in range(number_of_weeks) )-1)
            penalty = 2
            obj_int_vars.append(weekends_in_month[e])
            obj_int_coeffs.append(penalty)
            total_weekends[e] = model.NewIntVar(0, number_of_weeks+8,  'total_weekends[%i]' % e )    
            model.Add(total_weekends[e] == ansatte_number_of_weekends_worked_8weeks[ansatte[e]][0]+sum(weekenddays_worked[e, w] for w in range(number_of_weeks) ))
           # model.Add(total_weekends[e] == sum(weekenddays_worked[e, w] for w in range(number_of_weeks) ))
    
    #equalize number of weekends with work over 12 weeks (current+8 before current month)    
        min_fair_weekends_12weeks = model.NewIntVar(0, number_of_weeks+8, 'min_fair_shift__12weeks_%i' % e)
        max_fair_weekends_12weeks = model.NewIntVar(0, number_of_weeks+8, 'max_fair_shift__12weeks_%i' % e)
        model.AddMinEquality(min_fair_weekends_12weeks, [total_weekends[e] for e in range(num_employees)])        
        model.AddMaxEquality(max_fair_weekends_12weeks, [total_weekends[e] for e in range(num_employees)])
        weekends_diff_12weeks = model.NewIntVar(0, number_of_weeks+8, '')
        model.Add(weekends_diff_12weeks==max_fair_weekends_12weeks - min_fair_weekends_12weeks)
        penalty = 4
        obj_int_vars.append(weekends_diff_12weeks)
        obj_int_coeffs.append(penalty)
        
        

            
#fordel nattevagter via equality for alle nattevagter og fordel weekender
    sum_of_shifts = {}
    
    sum_of_equalfunctions = {}
    sum_of_equalshifts = {}
    for e in range(num_employees):
#nigt_shifts        
            sum_of_shifts[e] = model.NewIntVar(0, num_days, 'sum_of_shifts_%i' % e)
            model.Add(sum_of_shifts[e] == sum(work[e, s, d] for s in night_shifts for d in range(offset, num_days)))
#equalfunctions
            sum_of_equalfunctions[e] = model.NewIntVar(0, num_days, 'sum_of_equalfunctions_%i' % e)
            model.Add(sum_of_equalfunctions[e] == sum(work[e, s, d] for s in equalfunctions for d in range(offset, num_days)))
            
#night_shifts distrubuted equally across all nightshifts:
    if 2>1:
        sum_of_shifts_night = {}
        for e in range(num_employees):
            sum_of_shifts_night[e] = model.NewIntVar(0, num_days, 'sum_of_shifts_night_%i' % e)
            model.Add(sum_of_shifts_night[e] == sum(work[e, s, d] for s in night_shifts for d in range(offset, num_days)))
        
        for s in night_shifts:
            min_fair_shift_night = model.NewIntVar(0, num_days, 'min_fair_shift_night_%i' % s)
            max_fair_shift_night = model.NewIntVar(0, num_days, 'max_fair_shift_night_%i' % s)
            model.AddMinEquality(min_fair_shift_night, [sum_of_shifts_night[e] for e in range(num_employees)])
            model.AddMaxEquality(max_fair_shift_night, [sum_of_shifts[e] for e in range(num_employees)]) 
            nightshift_diff = model.NewIntVar(0, num_days, '')
            model.Add(nightshift_diff==max_fair_shift_night - min_fair_shift_night)
            penalty = 4
            obj_int_vars.append(nightshift_diff)
            obj_int_coeffs.append(penalty)

#equalfunctions
#    for s in equalfunctions:
#        min_fair_equalfunctions = model.NewIntVar(0, num_days, 'min_fair_equalfunctions_%i' % s)
#        max_fair_equalfunctions = model.NewIntVar(0, num_days, 'equalfunctions_%i' % s)
#        model.AddMinEquality(min_fair_equalfunctions, [sum_of_equalfunctions[e] for e in range(num_employees)])
#        model.AddMaxEquality(max_fair_equalfunctions, [sum_of_equalfunctions[e] for e in range(num_employees)]) 
#        equalfunctions_diff = model.NewIntVar(0, num_days, '')
#        model.Add(equalfunctions_diff==max_fair_equalfunctions - min_fair_equalfunctions)
#        penalty = 4
#        obj_int_vars.append(equalfunctions_diff)
#       obj_int_coeffs.append(penalty)

#equalfunctionsansatte
    sum_of_equalshifts = {}
    for s, ansatte_in_shift in equalfunctionsansatte:
        for e in ansatte_in_shift:
           
            sum_of_shifts[ansatte.index(e)] = model.NewIntVar(0, num_days, 'sum_of_shifts_%i' % ansatte.index(e))
            model.Add(sum_of_shifts[ansatte.index(e)] == sum(work[ansatte.index(e), s, d] for d in range(offset, num_days)))
        min_fair_equalfunctions = model.NewIntVar(0, num_days, 'min_fair_equalfunctions_%i' % s)
        max_fair_equalfunctions = model.NewIntVar(0, num_days, 'max_fair_equalfunctions_%i' % s)
        model.AddMinEquality(min_fair_equalfunctions, [sum_of_equalfunctions[ansatte.index(e)] for e in ansatte_in_shift])
        model.AddMaxEquality(max_fair_equalfunctions, [sum_of_equalfunctions[ansatte.index(e)] for e in ansatte_in_shift]) 
        equalfunctions_diff = model.NewIntVar(0, num_days, '')
        model.Add(equalfunctions_diff==max_fair_equalfunctions - min_fair_equalfunctions)
        penalty = 12
        obj_int_vars.append(equalfunctions_diff)
        obj_int_coeffs.append(penalty)



    
#limit workhours     
    month_end_hours = {
        e: model.NewIntVar(int(normtid*0.01* ansat_arbejdstid[e][1]/37), int(normtid*1.2 * ansat_arbejdstid[e][1]/37), 'employee_%s_month_end_hours' % e)
        for e in range(num_employees)
         
     }

    for e in range(num_employees):
        # this will be constrained because onth_end_hours is in domain [0, 24_800]
        tmp = sum( work[e, s, d] * int(funktioner_varighed[str(s)]) for s in range(1, num_shifts) for d in range (offset, num_days) )   
        model.Add(month_end_hours[e] == tmp)



 
# skills
    for e, s, k in skills:
        if k==0:
            for d in range (offset, num_days): #avoid conflicts by skillschange with month change
                if [e,s,d] not in fixed_assignments:#avoid fixed assignment conflicts
                #print ('e',e,'s',s,'k',k, 'd', d) 
                    model.Add(work[e, s, d] == 0)


#enforce fkjt 2 efter fkt1
#    for e in range(num_employees):
#              for d in range(num_days-1):
#                  model.Add(work[e, 2, d + 1] == 1 ).OnlyEnforceIf(work[e, 1, d])   
#enforce day shifts
#    for e in range(num_employees):
#        for w in range(num_weeks):
#            model.Add(work[e, 2,  w*7+6] == 1 ).OnlyEnforceIf(work[e, 1, w*7+5])
            


#fri hele weekenden                 
    if 1>2:
        for e in range(num_employees):
            for w in range(num_weeks):
                desired = [
                            work[e, 0, 7*w+6].Not(), 
                            work[e, 0, 7*w+5] 
                                ]
                desired_var = model.NewBoolVar(
                        'whole weekend free (employee=%i, w=%i)' % (e, w))
                desired.append(desired_var)
                model.AddBoolOr(desired)
                obj_bool_vars.append(desired_var)
                obj_bool_coeffs.append(2)               


    #     for w in range(num_weeks):
#             model.Add(work[e, 0, 7*w+6] == 1 ).OnlyEnforceIf(work[e, 0, 7*w+5])            

# favor desired_day_transistions
    for first_shift, first_weekday, second_shift, second_weekday, cost in desired_day_transistions:
        for e in range(num_employees):
            for w in range(num_weeks):
                desired = [
                                work[e, first_shift, (7*w + first_weekday)].Not(), 
                                work[e, second_shift, (7*w + second_weekday)] 
                ]
                #model.AddBoolOr(desired)
                if cost == 0:
                    #model.AddBoolOr(desired)
                    model.Add(work[e, second_shift, (7*w + second_weekday)] == 1 ).OnlyEnforceIf(work[e, first_shift, (7*w + first_weekday)])
                else:
                    desired_var = model.NewBoolVar(
                            'desired_day_transition (employee=%i, day=%i, shift=%i w=%i)' % (e, 7*w + first_weekday, first_shift, cost))
                    desired.append(desired_var)
                    model.AddBoolOr(desired)
                    obj_bool_vars.append(desired_var)
                    obj_bool_coeffs.append(-1*cost)               

# rewarded - desired_shift_transitions
    #for previous_shift, next_shift, cost in rewarded_transitions:
    for first_shift, first_day, next_shift, next_day, cost in desired_shift_transitions:
        for e in range(num_employees):
            for d in range((num_days - next_day)):
                transition_wanted = [
                    work[e, first_shift, d+first_day].Not(), work[e, next_shift,
                                                           d + next_day]
                ]
                
                if cost == 0:
                    #model.AddBoolOr(transition_wanted)
                    model.Add(work[e, next_shift, d + next_day] == 1 ).OnlyEnforceIf(work[e, first_shift, d+first_day])
                #if cost == 4:    
                    #model.Add(work[e, next_shift, d + next_day] == 1 ).OnlyEnforceIf(work[e, first_shift, d+first_day])
                    
                    
                else:
                    trans_wanted_var = model.NewBoolVar(
                        'desired_shift_transition (employee=%i, day=%i, firstshift=%i)' % (e, d, first_shift))
                    transition_wanted.append(trans_wanted_var)
                    model.AddBoolOr(transition_wanted)
                    obj_bool_vars.append(trans_wanted_var)
                    obj_bool_coeffs.append(-1*cost)
    
#only one weekend per month    
    if 1>2:
       
            cost=4
        
            for week in range(num_weeks):
                print('week',week)
                print('break')
                for w in range (week,num_weeks):
                    if w>week:
                        print ('w',w)
                        for e in range(num_employees):
                            #lør+søn
                            for s in funktioner_vagt_index_sleep: 
                            #+funktioner_vagt_dag_index_sleep:
                                #print ('w',w)
                                trans_wanted_sat=[work[e, s, week*7+5].Not(), work[e, 0,  w*7+5]]
                                trans_wanted_sat2=[work[e, s, week*7+5].Not(), work[e, 0,  w*7+6]]
                                trans_wanted_sun=[work[e, s, week*7+6].Not(), work[e, 0, w*7+6]]
                                trans_wanted_sun2=[work[e, s, week*7+6].Not(), work[e, 0, w*7+5]]
                                trans_wanted_sat_var=model.NewBoolVar('desired_shift_transition_week__weekend_sat (employee=%i, week=%i, weekend=%i firstshift=%i)' % (e, week, w, s))
                                trans_wanted_sat2_var=model.NewBoolVar('desired_shift_transition_week__weekend_sat2 (employee=%i, week=%i, weekend=%i firstshift=%i)' % (e, week, w, s))
                                trans_wanted_sun_var=model.NewBoolVar('desired_shift_transition_week__weekend_sun (employee=%i, week=%i, weekend=%i firstshift=%i)' % (e, week, w, s))
                                trans_wanted_sun2_var=model.NewBoolVar('desired_shift_transition_week__weekend_sun2 (employee=%i, week=%i, weekend=%i firstshift=%i)' % (e, week, w, s))    
                                trans_wanted_sat.append(trans_wanted_sat_var)
                                trans_wanted_sat2.append(trans_wanted_sat2_var)
                                trans_wanted_sun.append(trans_wanted_sun_var)
                                trans_wanted_sun2.append(trans_wanted_sun2_var)
                                model.AddBoolOr(trans_wanted_sat)
                                model.AddBoolOr(trans_wanted_sat2)
                                model.AddBoolOr(trans_wanted_sun)
                                model.AddBoolOr(trans_wanted_sun2)
                                obj_bool_vars.append(trans_wanted_sat_var)
                                obj_bool_coeffs.append(cost)
                                obj_bool_vars.append(trans_wanted_sat2_var)
                                obj_bool_coeffs.append(cost)
                                obj_bool_vars.append(trans_wanted_sun_var)
                                obj_bool_coeffs.append(cost)
                                obj_bool_vars.append(trans_wanted_sun2_var)
                                obj_bool_coeffs.append(cost)
                            #fredag nat    
                            for s in night_shifts:
                                cost=4
                                trans_wanted_fri=[work[e, s, week*7+4].Not(), work[e, 0,  w*7+4]]
                                trans_wanted_fri_var=model.NewBoolVar('desired_shift_transition_week_fri (employee=%i, week=%i, weekend=%i firstshift=%i)' % (e,week, w, s))
                                trans_wanted_fri.append(trans_wanted_fri_var)     
                                model.AddBoolOr(trans_wanted_fri)
                                obj_bool_vars.append(trans_wanted_fri_var)
                                obj_bool_coeffs.append(cost)
                                trans_wanted_fri2=[work[e, s, week*7+4].Not(), work[e, 0,  w*7+5]]
                                trans_wanted_fri2_var=model.NewBoolVar('desired_shift_transition_week_fri2 (employee=%i, week=%i, weekend=%i firstshift=%i)' % (e,week, w, s))
                                trans_wanted_fri2.append(trans_wanted_fri2_var)     
                                model.AddBoolOr(trans_wanted_fri2)
                                obj_bool_vars.append(trans_wanted_fri2_var)
                                obj_bool_coeffs.append(cost)
                                trans_wanted_fri3=[work[e, s, week*7+4].Not(), work[e, 0,  w*7+6]]
                                trans_wanted_fri3_var=model.NewBoolVar('desired_shift_transition_week_fri3 (employee=%i, week=%i, weekend=%i firstshift=%i)' % (e,week, w, s))
                                trans_wanted_fri3.append(trans_wanted_fri3_var)     
                                model.AddBoolOr(trans_wanted_fri3)
                                obj_bool_vars.append(trans_wanted_fri3_var)
                                obj_bool_coeffs.append(cost)
                    

#avoid 2 consecutive weekends
    if 1<2:
        for e in range(num_employees):
            for w in range(num_weeks-2):
                #sat+sun
                for s in funktioner_vagt_index_sleep: 
                #+funktioner_vagt_dag_index_sleep:
                    cost=8
                    trans_wanted_sat=[work[e, s, w*7+5].Not(), work[e, 0,  w*7+12]]
                    trans_wanted_sat_var=model.NewBoolVar('desired_shift_transition_week_sat (employee=%i, week=%i, firstshift=%i)' % (e, w, s))
                    trans_wanted_sat.append(trans_wanted_sat_var)     
                    model.AddBoolOr(trans_wanted_sat)
                    obj_bool_vars.append(trans_wanted_sat_var)
                    obj_bool_coeffs.append(cost)
                    trans_wanted_sun=[work[e, s, w*7+6].Not(), work[e, 0, w*7+13]]
                    trans_wanted_sun_var=model.NewBoolVar('desired_shift_transition_week_sun (employee=%i, week=%i, firstshift=%i)' % (e, w, s))
                    trans_wanted_sun.append(trans_wanted_sun_var)     
                    model.AddBoolOr(trans_wanted_sun)
                    obj_bool_vars.append(trans_wanted_sun_var)
                    obj_bool_coeffs.append(cost)
        
                    trans_wanted_sat2=[work[e, s, w*7+5].Not(), work[e, 0,  w*7+13]]
                    trans_wanted_sat2_var=model.NewBoolVar('desired_shift_transition_week_sat2 (employee=%i, week=%i, firstshift=%i)' % (e, w, s))
                    trans_wanted_sat2.append(trans_wanted_sat2_var)     
                    model.AddBoolOr(trans_wanted_sat2)
                    obj_bool_vars.append(trans_wanted_sat2_var)
                    obj_bool_coeffs.append(cost)
                    trans_wanted_sun2=[work[e, s, w*7+6].Not(), work[e, 0, w*7+12]]
                    trans_wanted_sun2_var=model.NewBoolVar('desired_shift_transition_week_sun2 (employee=%i, week=%i, firstshift=%i)' % (e, w, s))
                    trans_wanted_sun2.append(trans_wanted_sun2_var)     
                    model.AddBoolOr(trans_wanted_sun2)
                    obj_bool_vars.append(trans_wanted_sun2_var)
                    obj_bool_coeffs.append(cost)
                #friday lateshifts
                for s in night_shifts:
                    cost=8
                    trans_wanted_fri=[work[e, s, w*7+4].Not(), work[e, 0,  w*7+11]]
                    trans_wanted_fri_var=model.NewBoolVar('desired_shift_transition_week_fri (employee=%i, week=%i, firstshift=%i)' % (e, w, s))
                    trans_wanted_fri.append(trans_wanted_fri_var)     
                    model.AddBoolOr(trans_wanted_fri)
                    obj_bool_vars.append(trans_wanted_fri_var)
                    obj_bool_coeffs.append(cost)
                    trans_wanted_fri2=[work[e, s, w*7+4].Not(), work[e, 0, w*7+12]]
                    trans_wanted_fri2_var=model.NewBoolVar('desired_shift_transition_week_fri2 (employee=%i, week=%i, firstshift=%i)' % (e, w, s))
                    trans_wanted_fri2.append(trans_wanted_fri2_var)     
                    model.AddBoolOr(trans_wanted_fri2)
                    obj_bool_vars.append(trans_wanted_fri2_var)
                    obj_bool_coeffs.append(cost)
                    
                    trans_wanted_fri3=[work[e, s, w*7+4].Not(), work[e, 0, w*7+13]]
                    trans_wanted_fri3_var=model.NewBoolVar('desired_shift_transition_week_fri3 (employee=%i, week=%i, firstshift=%i)' % (e, w, s))
                    trans_wanted_fri3.append(trans_wanted_fri3_var)     
                    model.AddBoolOr(trans_wanted_fri3)
                    obj_bool_vars.append(trans_wanted_fri3_var)
                    obj_bool_coeffs.append(cost)
        
                    
# Shift constraints
    for ct in shift_constraints:
        shift, hard_min, soft_min, min_cost, soft_max, hard_max, max_cost = ct
        for e in range(num_employees):
            works = [work[e, shift, d] for d in range(num_days)]
            variables, coeffs = add_soft_sequence_constraint(
                model, works, hard_min, soft_min, min_cost, soft_max, hard_max,
                max_cost, 'shift_constraint(employee %i, shift %i)' % (e,
                                                                       shift))
            obj_bool_vars.extend(variables)
            obj_bool_coeffs.extend(coeffs)

 #evenly distribution
#    max_days_per_e = num_days #((num_days//7)*2)-5 
#    min_days_per_e = 0 #5 #num_days-((num_days//7)*2)    #num_days
    
    #print ('num_days', num_days,'min', min_days_per_e, 'max', max_days_per_e, 'n e', num_employees )
#    for e in range(num_employees):
#        num_days_worked = sum(work[e, s, d] for s in range(1, num_shifts) for d in range(num_days))#- shift 0 off
#        model.Add(min_days_per_e <= num_days_worked)
#        model.Add(num_days_worked <= max_days_per_e)
#        #print (num_days_worked)





# Objective
    model.Minimize(
        sum(obj_bool_vars[i] * obj_bool_coeffs[i]
            for i in range(len(obj_bool_vars)))
        + sum(obj_int_vars[i] * obj_int_coeffs[i]
              for i in range(len(obj_int_vars))))
       


    if output_proto:
        print('Writing proto to %s' % output_proto)
        with open(output_proto, 'w') as text_file:
            text_file.write(str(model))

    # Solve the model.
    
    solver = cp_model.CpSolver()
    solver.parameters.max_time_in_seconds = calc_time
    solver.parameters.num_search_workers = 8 
    if params:
        text_format.Merge(params, solver.parameters)
    solution_printer = cp_model.ObjectiveSolutionPrinter()
    status = solver.SolveWithSolutionCallback(model, solution_printer)

    # Print solution.
    if status == cp_model.OPTIMAL or status == cp_model.FEASIBLE:
        if 2>1:
            print()
            a=''
            for e in range(num_employees):
                for w in range(number_of_weeks):
                    #b=solver.Value(weekends[e, w])
                    y=solver.Value(weekenddays_worked_int[e, w])
                    p=solver.BooleanValue(weekenddays_worked[e, w])
                    a+='e%i_w%i_weekdays%i_bool%i' %(e,w,y,p) + '\n'
                    
                z=solver.Value(total_weekends[e]) 
                a+='empl %s total weekends %i' %(ansatte[e],z)+ '\n'
            print(a)
        print()
        header = '  '
        for w in range(num_weeks):
            header += '  M    T    W    T    F    S    S   '
        print(header)
        for e in range(num_employees_dummy):
            schedule = ''
            for d in range(offset, num_days):
                for s in range(num_shifts):
                    if solver.BooleanValue(work[e, s, d]):
                        schedule += shifts[s] + ' '
            ny_e =  "{0:0>2}".format(e)
          #  print(employees[e])  
           # print('worker:',ny_e, '%s' % (schedule))
            print(employees_dummy[e]+':', '%s' % (schedule))

        print()
        allskema=""
        for d in range(num_days):
            day = d
            date = dates[d]
            ugedag = ugedag_dag[d]
            skema = '{'+ f'day : {day} '+f'-{ugedag}'+f' ({date}) :'
            fri = ', fri : "'
            for s in range(num_shifts): 
                for e in range(num_employees_dummy):
                    if solver.BooleanValue(work[e, s, d]):
                        if shifts[s] != "Sove":
                            skema += shifts[s] + ':' +  employees_dummy[e]+','
                        if  shifts[s] == "Sove":
                            fri += employees_dummy[e] + ','
            skema +=  '"'+fri[:-1]+'"'+'}'       
            allskema += skema + '\n'           
            print(skema)
        print()
        data={}
        for s in range(num_shifts):
            data[shifts[s]]=[]
            for d in range(num_days):
                ansatte =''
                for e in range(num_employees_dummy):
                    if solver.BooleanValue(work[e, s, d]):
                        ansatte += employees_dummy[e]+" "
                        #data[shifts[s]].append({dates[d]:employees[e]})
                data[shifts[s]].append({dates[d]:ansatte.rstrip()})
        import json
        with open('brain.json', 'w') as outfile:
            json.dump(data, outfile)                
                        
        
        print(data)
                    
        
        
            
            
        outF = open("OutFile.txt", "w")
        outF.write(allskema)
        outF.close()
            
        print()
        print('Penalties:')
        for i, var in enumerate(obj_bool_vars):
            if solver.BooleanValue(var):
                penalty = obj_bool_coeffs[i]
                if penalty > 0:
                    print('  %s violated, penalty=%i' % (var.Name(), penalty))
                else:
                    print('  %s fulfilled, gain=%i' % (var.Name(), -penalty))

        for i, var in enumerate(obj_int_vars):
            if solver.Value(var) > 0:
                print('  %s violated by %i, linear penalty=%i' %
                      (var.Name(), solver.Value(var), obj_int_coeffs[i]))
                
    print()
    
    print(solver.ResponseStats())
    #print(obj_bool_vars) 

def main(args):
    """Main."""
    solve_shift_scheduling(args.params, args.output_proto)


if __name__ == '__main__':
    main(PARSER.parse_args())
