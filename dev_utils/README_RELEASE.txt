
Update release number in 
	poppy/_version.py
	docsource/installation file names and links (TBD: update programmatically?)
	make-data-sdist.sh


Make the sdists
	python setup.py sdist

Tag the release:
	svn copy svn+ssh://tib.stsci.edu/home/mperrin/Repository/svnbase/jwst/webbpsf/trunk svn+ssh://tib.stsci.edu/home/mperrin/Repository/svnbase/jwst/webbpsf/tags/release-0.2.6

