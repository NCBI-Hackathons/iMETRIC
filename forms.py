from flask_wtf import Form
import wtforms
from wtforms.fields import *
from wtforms.widgets import *

class SubmissionForm(Form):
    protein_sequence = TextAreaField('Protein Sequence',
                            validators=[validators.input_required()])
    submit = SubmitField("Send")
