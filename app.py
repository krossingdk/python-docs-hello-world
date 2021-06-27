#from flask import Flask
#app = Flask(__name__)

#@app.route("/")
#def hello():
#    return "Hello, World!"


import os
import sys
import database
import pickle
from flask import Flask, request, render_template, redirect, url_for
import json



project_root = os.path.dirname(os.path.realpath('__file__'))
template_path = os.path.join(project_root, 'app/templates')
static_path = os.path.join(project_root, 'app/static')
app = Flask(__name__, template_folder=template_path, static_folder=static_path)


@app.route('/')
def index():
    #exec(open('schedule.py').read())
    return 'On the way to automatic scheduling'


@app.route('/vagtplan')
def vagtplan():
    import os
    os.system("python schedule.py")
    file = open('OutFile.txt')
    txt= file.read()
    
    return txt 

@app.route('/basic')
def basic():
    import basic
    return 'Hello running basic task...'

@app.route('/ansatte')
def ansatte():
   
    return database.ansatte()

@app.route('/titel')
def titel():
    return database.titel()

@app.route('/funktioner')
def funktioner():
    return database.funktioner()

@app.route('/json', methods=['POST'])
def json():
    
    import json
    data={}
    with open('brain.json', 'w') as outfile:
            json.dump(data, outfile)
    request_data = request.get_json()
    with open('data.json', 'w') as file:
        json.dump(request_data ,file)
    import os
    os.system("python brain2.py")
    with open('brain.json') as json_file:
        data = json.load(json_file)
    #return request_data
    return data

@app.route('/reply', methods=['GET', 'POST'])
def reply():
    skills = request.args.get('skills')
    requests = request.args.get('requests')
    cover_demands = request.args.get('cover_demands')
    weekly_sum_constraints = request.args.get('weekly_sum_constraints')
    funktioner_id_all = request.args.get('funktioner_id_all')
    funktioner_funktion_all = request.args.get('funktioner_funktion_all')
    funktioner_dag_index_sleep = request.args.get('funktioner_dag_index_sleep')
    funktioner_vagt_index_sleep = request.args.get('funktioner_vagt_index_sleep')
    ansatte = request.args.get('ansatte')
    offset = request.args.get('offset')
    period_length = request.args.get('period_length')
    number_of_weeks = request.args.get('number_of_weeks')
    fixed_assignments = request.args.get('fixed_assignments')
    
    with open('data.pkl', 'wb') as file:
        pickle.dump(skills ,file)
        pickle.dump(requests,file)
        pickle.dump(cover_demands,file)
        pickle.dump(weekly_sum_constraints,file)
        pickle.dump(funktioner_id_all,file)
        pickle.dump(funktioner_funktion_all,file)
        pickle.dump(funktioner_dag_index_sleep,file)
        pickle.dump(funktioner_vagt_index_sleep,file)
        pickle.dump(ansatte,file)
        pickle.dump(offset,file)
        pickle.dump(period_length,file)
        pickle.dump(number_of_weeks,file)
        pickle.dump(fixed_assignments,file)
        pickle.dump(d,file)
    
    
    d = dict();  
    d['skills'] = skills
    d['requests']   = requests
    d['cover_demands']   = cover_demands
    d['weekly_sum_constraints']   = weekly_sum_constraints
    d['funktioner_id_all'] = funktioner_id_all
    d['funktioner_funktion_all'] = funktioner_funktion_all
    d['funktioner_dag_index_sleep'] = funktioner_dag_index_sleep
    d['funktioner_vagt_index_sleep'] = funktioner_vagt_index_sleep
    d['ansatte'] = ansatte
    d['offset'] = offset
    d['period_length'] = period_length
    d['number_of_weeks'] = number_of_weeks
    d['fixed_assignments'] = fixed_assignments 
    
    
    
    
    import os
    os.system("python brain.py")    

   
    
    return d 

    
application = app

if __name__ == "__main__":
    app.run()
