
To upload a new version to PYPI (as two separate packages):

	# actually this would work better if we did
	
	\cp setup-poppy.py setup.py
	python setup.py sdist upload

	\cp setup-webbpsf.py setup.py
	python setup-webbpsf.py sdist upload

To create a unified tar file containing both codes:

	\cp setup-both.py setup.py
	python setup-webbpsf.py sdist upload

To create a distribution file for the data: 

	
