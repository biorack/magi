
import requests as re

"""
Integrate MAGI outputs to easily order DNA through JGI BOOST
https://boost.jgi.doe.gov/boost.html

This is a work in progress and not ready to be used yet.

Our plan is make a function here that takes these inputs
* a set of sequences
* boost-token, and 
* boost-parameters
This will create an entry in BOOST that can be sent to 
the JGI DNA synthesis group in order to build it.
"""

# STEP 1.
# sign into boost with jgi SSO

# STEP 2.
# get key from https://boost.jgi.doe.gov/boost.html#user
# copy jwt into my_key below
my_key = 'my_key'

# STEP 3.
# Define workflow as json text
#This is an example workflow from 
# https://boost.jgi.doe.gov/boost.html#manual
# BOOST API
# Submission of Jobs
# Submission of a Workflow including all BOOST functionalities.
with open('/Users/bpb/Downloads/example_boost_workflow_submission.json','r') as fid:
    job_description = fid.read()
    
# STEP 4.
#Submit this job to boost along with your key
r = re.post('https://boost.jgi.doe.gov/rest/jobs/submit',
            data = job_description,
            headers={'content-type': 'application/json',
                     'authorization':my_key})

# STEP 5.
# Get job UUID in the response to POST request
print(r.text)

# STEP 6. Repeat asking API until r.json()['job']['job-status']=='FINISHED'
job_uuid = 'put uuid here'
r = re.get('https://boost.jgi.doe.gov/rest/jobs/%s'%job_uuid,
            headers={'content-type': 'application/json',
                     'authorization':my_key})

# STEP 7.
# see if its running
print(r.text)
