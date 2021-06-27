
from flask import Flask
import json

app = Flask(__name__)


@app.route('/')
def index():
    return 'On the way to automatic scheduling'



