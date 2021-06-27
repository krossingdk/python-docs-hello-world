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
