from flask import Flask
app = Flask(__name__)
import os
import sys
import database
import pickle
import json

@app.route("/")
def hello():
    return "Hello, Kasper Rossing!"

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
