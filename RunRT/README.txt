#To create a stand-alone python app for mac run pyinstaller

pyinstaller -w RunRT.py

cp *.pbm *.dat runrtdoc.txt rtdoc.txt dist/RunRT
cp -r RUNS dist/RunRT
