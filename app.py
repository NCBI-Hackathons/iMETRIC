#System Imports 
import sys, os
import json
import static
import time
import random
from shutil import copyfile
import operator
import urllib2
import itertools
import subprocess
import math 
# from filechunkio import FileChunkIO 
from celery import Celery
from collections import defaultdict, OrderedDict
import collections
import requests
#Flask Imports
from werkzeug import secure_filename
from flask import Flask, Blueprint, make_response, render_template, render_template_string, request, session, flash, redirect, url_for, jsonify, get_flashed_messages, send_from_directory
from flask.ext.bcrypt import Bcrypt
from flask.ext.login import LoginManager, UserMixin, current_user, login_user, logout_user, login_required
from flask.ext.mail import Mail, Message
from flask.ext.script import Manager
from flask.ext.migrate import Migrate, MigrateCommand
from flask_bootstrap import Bootstrap
from flask_bootstrap import __version__ as FLASK_BOOTSTRAP_VERSION
from flask_nav import Nav 
from flask_nav.elements import Navbar, View, Subgroup, Link, Text, Separator
from flask_sqlalchemy import SQLAlchemy
from markupsafe import escape
import wtforms
from flask_wtf import Form
import random
import jinja2 
from flask_restful import reqparse, abort, Api, Resource

from sqlalchemy import create_engine
from sqlalchemy import Table, Column, Integer, String, MetaData, ForeignKey, Boolean
from sqlalchemy.dialects import postgresql
from sqlalchemy.dialects.postgresql import JSON, JSONB, ARRAY, BIT, VARCHAR, INTEGER, FLOAT, NUMERIC, OID, REAL, TEXT, TIME, TIMESTAMP, TSRANGE, UUID, NUMRANGE, DATERANGE
from sqlalchemy.sql import select
from sqlalchemy.orm import sessionmaker, scoped_session



app = Flask(__name__)
api = Api(app)

# Get appropriate port from Heroku (or where ever we host)
try:
 port = int(os.environ['PORT'])
except:
 port = 5001
print 'will run on port ' + str(port)


@app.route('/')
def index():
    return render_template("index.html")




# DATA-GATHERING FUNCTIONS 

def post_to_iedb_mhci(protein_sequence, method='smm', length='9', allele='HLA-A*01:01'): 
	data = {
	'sequence_text': protein_sequence, 
	'length': length,
	'method': method, 
	'allele': allele,
	}
	url = 'http://tools-api.iedb.org/tools_api/mhci/'
	response = requests.post(url, data=data) 
	if response.ok: 
		return response.text 
	else: 
		return 'Something went wrong'



# API RESOURCES BELOW 

class ProteinQuery(Resource): 
    def get(self, protein_sequence): 
			response = post_to_iedb_mhci(protein_sequence)
			return {protein_sequence: response}

api.add_resource(ProteinQuery, '/query/<string:protein_sequence>')


if __name__ == '__main__':
    app.run(port=port, debug=True)
