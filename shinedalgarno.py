import json
import os
import subprocess
import re
import progressbar


def calculate_values_for_sequences(sequences: list):
  print(os.getcwd())
  resp = subprocess.Popen(r"..\tools\ViennaRNA\RNAcofold.exe ../data/shinedalgarno.txt", shell=True,
                          stdout=subprocess.PIPE)
  resp = resp.stdout.read().decode().splitlines()
  lines = list(filter(lambda l: l.strip().endswith(')'), resp))
  values = list(map(float, re.findall(r'\s\(([^\)]+)\)', "\n".join(l[6:] for l in lines))))

  return [dict(seq=var, value=val) for var, val in zip(sequences, values)]


def calculate_from_file_and_save_to_json():
  with open('../data/shinedalgarno.txt') as fp:
    sequences = [ln.strip() for ln in fp.readlines()]

  print('calculating')
  shinedal = calculate_values_for_sequences(sequences)
  print(shinedal)
  # create folds.json file
  with open('../data/shinedalgarno.json', 'w', encoding='utf8') as fp:
    json.dump(shinedal, fp, indent=2, ensure_ascii=False)

if __name__ == '__main__':
  calculate_from_file_and_save_to_json()
