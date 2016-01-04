import sys, os
from flask import Flask
from flask_restful import reqparse, abort, Api, Resource
from flask import render_template

app = Flask(__name__)
api = Api(app)

# Get appropriate port from Heroku (or where ever we host)
try:
 port = int(os.environ['PORT'])
except:
 port = 5001

print 'running on port ' + str(port)

@app.route('/')
def index():
    return render_template("index.html")

if __name__ == '__main__':
    app.run(port=port, debug=True)
