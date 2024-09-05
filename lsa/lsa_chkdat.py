#!/usr/bin/env python
import argparse
import numpy as np


def main():
  # check_data datafile.txt
  parser = argparse.ArgumentParser(description="Auxillary tool to new LSA package for checking data format")
  parser.add_argument("dataFile", metavar= "dataFile", type=argparse.FileType('r'), help="the data file")
  parser.add_argument("repNum", metavar= "repNum", type=int, help="replicates number")
  parser.add_argument("spotNum", metavar= "spotNum", type=int, help="timepoints number")

  arg_namespace = parser.parse_args()

  dataFile = vars(arg_namespace)['dataFile']
  repNum = vars(arg_namespace)['repNum']
  spotNum = vars(arg_namespace)['spotNum']
  colNum = repNum*spotNum

  index = 0
  for row in dataFile:
    if row[0]== '#':
      index += 1
      continue
    else:
      cells= row.strip('\n').split('\t')
      try:
        assert len(cells) == colNum + 1
      except AssertionError:
        print("Format error:", index, "-th row has",len(cells),"cells instead of expected number:", colNum+1)
        quit(0)
      for i in range(1, colNum+1):
        try:
          float(cells[i])
        except:
          try:
            assert cells[i] in ['na','NA','']
          except AssertionError:
            print()
      index += 1
  dataFile.seek(0)
  try:
    np.genfromtxt( \
        dataFile, comments='#', delimiter='\t', missing_values=['na','','NA'], \
        filling_values=np.nan, usecols=list(range(1,spotNum*repNum+1)) )
  except ValueError:
    print("Numpy cann't read", str(ValueError))

if __name__=="__main__":
  main()
