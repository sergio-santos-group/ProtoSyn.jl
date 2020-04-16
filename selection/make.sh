antlr4='java -jar /usr/local/lib/antlr-4.7.1-complete.jar'

# JAVA
$antlr4 Selection.g4  -no-listener -o java && javac java/*java


# PYTHON
#sudo pip install antlr4-python2-runtime
#python main.py ../tmp2.pbt



#$antlr4 -Dlanguage=Python3 Selection.g4 -visitor -no-listener -o src
#javac Selection*.java

