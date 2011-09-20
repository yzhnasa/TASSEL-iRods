@echo off

set TOP=.
set LIB_JARS=%TOP%\lib

set CP=

set CP=.\dist\sTASSEL.jar
for %%i in (%LIB_JARS%\*.jar) do call ".\cp.bat" %%i
echo %CP%

java -classpath "%CP%" -Xms512m -Xmx1536m net.maizegenetics.tassel.TASSELMainApp %*
