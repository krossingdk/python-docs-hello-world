

import os
import sys
import database
import pickle
from flask import Flask, request, render_template, redirect, url_for
import json

app = Flask(__name__)


@app.route('/')
def index():
    return 'On the way to automatic scheduling'



