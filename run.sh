#!/bin/bash

set -e

here=`pwd`
KATIEpath=/root/physics/katieWithPDFxTMD
create_lib=$KATIEpath/create_lib.py
create_hef=$KATIEpath/create_hef.py
build=$KATIEpath/build

task=""
sourceFile=""
executable=""
files=""
topic=""

function info {
  echo "=============================================================================="
  echo " Examples: (abbreviating '$0' as 'run.sh')"
  echo ""
  echo " $ run.sh clean"
  echo " $ run.sh lib"
  echo " $ run.sh prepare <filename> <dirname>"
  echo " $ run.sh compile <sourcefile>"
  echo " $ run.sh compile,run <sourcefile>"
  echo " $ run.sh compile,run <sourcefile> <datafile>"
  echo " $ run.sh merge raw1.dat raw2.dat raw3.dat"
  echo " $ run.sh merge raw*"
  echo " $ run.sh lhef raw1.dat raw2.dat raw3.dat"
  echo " $ run.sh lhef raw*"
  echo " $ run.sh help compile"
  echo " $ run.sh katamp"
  echo " "
  echo " You can also give arguments with explicit keywords."
  echo " Then the arguments do not need to follow a specific order."
  echo " "
  echo " $ run.sh task=clean"
  echo " $ run.sh task=lib"
  echo " $ run.sh task=prepare input=<filename> dir=<dirname>"
  echo " $ run.sh task=compile source=<sourcefile>"
  echo " $ run.sh task=compile,run source=<sourcefile>"
  echo " $ run.sh task=compile,run source=<sourcefile> files=<datafile>"
  echo " $ run.sh task=merge files=raw1.dat,raw2.dat,raw3.dat"
  echo " $ run.sh task=help topic=compile"
  echo " $ run.sh task=katamp"
  echo "=============================================================================="
}

if [ "$#" -eq 0 ]; then
  info
  exit 0
fi

for var in "$@"
do
    lhs=`echo $var | sed s~=.*$~~`
    rhs=`echo $var | sed s~^$lhs=~~`
    case $lhs in
    "task")
        task=$rhs
        ;;
    "seed")
        seed=$rhs
        ;;
    "Nev")
        Nev=$rhs
        ;;
    "input")
        input=$rhs
        ;;
    "dir")
        targetDir=$rhs
        ;;
    "source")
        sourceFile=`echo $rhs | sed -e's/,/ /g'`
        ;;
    "exec")
        executable=$rhs
        ;;
    "files")
        files=`echo $rhs | sed -e's/,/ /g'`
        ;;
    "topic")
        topic=$rhs
        ;;
    esac
done


if [ "$task" = "" ]; then
    task=$1
    taskList=(`echo "$task" | sed 's/,/ /g'`)
    case ${taskList[0]} in
    "prepare")
        input=$2
        targetDir=$3
        ;;
    "run")
        executable=$2
        Nev=$3
        seed=$4
        targetDir=$5
        ;;
    "compile")
        case ${taskList[1]} in
        "run")
            sourceFile=$2
            if [ "$#" -gt 2 ]; then
              files="${@:3}"
            fi
            ;;
        *)
            #sourceFile=$2
            sourceFile="${@:2}"
            ;;
        esac
        ;;
    "merge"|"lhef")
        files="${@:2}"
        ;;
    "help")
        topic=$2
        ;;
    esac
fi


taskList=(`echo "$task" | sed 's/,/ /g'`)
case ${taskList[0]} in
"clean")
    mkdir -p "$build"
    rm -f "$build"/*
    ;;
"lib")
    case ${taskList[1]} in
    "twostep")
        echo 'Writing source files...'
        python $create_lib lib print | grep executing \
               | sed -e's/executing: //' > $build/compile_all.sh
        echo 'Compiling source files...'
        bash $build/compile_all.sh
        ;;
    *)
        python $create_lib lib
        ;;
    esac
    ;;
"compile")
    case ${taskList[1]} in
    "run")
        python $create_lib compile $sourceFile
        sourceList=($sourceFile)
        ./`echo ${sourceList[-1]} | sed -e's/\.f[0-9][0-9]$/.out/'` $files
        ;;
    *)
        python $create_lib compile $sourceFile
        ;;
    esac
    ;;
"prepare")
    if [[ "$targetDir" != "./" && "$targetDir" != "." ]]; then
        mkdir $targetDir
    fi
    touch $targetDir/extra_cuts.h90
    touch $targetDir/extra_weights.h90
    python $create_hef "input="$input "dir="$targetDir
    cp $input $targetDir/input
    chmod 700 $targetDir/optimize.sh
    #chmod 700 $targetDir/recompile.sh
    #chmod 700 $targetDir/create_eventfile.sh
    ;;
"run")
    mkdir $targetDir
    if [ "$executable" = "" ]; then
        executable=main.out
    fi
    cp $executable $targetDir/$targetDir$executable
    cd $targetDir
        ./$targetDir$executable "seed="$seed "Nev="$Nev "dir=./" > output
        rm $targetDir$executable
    cd $here
    ;;
"merge"|"lhef")
    $build/merge_raw.out $task $files
    ;;
"help")
    if [ "$topic" = "" ]; then
      info
    else
      python $create_lib help $topic
    fi
    ;;
"katamp")
    python $create_lib $task
    ;;
"source")
    python $create_lib $task
    ;;
*)
    echo "ERROR in $0: task $task is not defined"
    exit 1
esac
