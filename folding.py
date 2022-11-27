import collections
import json
import subprocess
import re
import progressbar


def variant_to_frames(variant: str, frame_length=30):
  return [variant[i:i + frame_length] for i in range(len(variant) - frame_length + 1)]


def calculate_values_for_variants(variants: list):
  with open(r'temp_variants.txt', 'w') as fp:
    print(f'writing {len(variants)} variants to temp_variants.txt')
    fp.writelines([v+'\n' for v in variants])

  resp = subprocess.Popen(r".\ViennaRNA\RNAfold.exe < temp_variants.txt", shell=True, stdout=subprocess.PIPE)

  value_pattern = r'\(\s(-?\d+.\d+)\)'
  values = map(float, re.findall(value_pattern, resp.stdout.read().decode()))

  return [dict(variant=var, value=val) for var, val in zip(variants, values)]


def folding_for_list(variants):
  print(f'{len(variants)} sequences found')

  variants_folds = []

  for variant in progressbar.progressbar(variants):
    folds = calculate_values_for_variants([variant, *variant_to_frames(variant)])
    fold_main = folds[0]
    folds_local_values = [f['value'] for f in folds[1:]]

    variants_folds.append(dict(
      variant={k: v for k, v in folds[0].items() if k != 'frame'},
      localValues=folds_local_values
    ))

  return variants_folds


def folding_for_csv(csv_file: str):
  with open(csv_file) as fp:
    variants = list(filter(bool, (ln.strip() for ln in fp.readlines()[1:])))  # read as txt file
  return folding_for_list(variants)


if __name__ == '__main__':
  pass

  # folding_for_csv('../data/Variants.csv')
  #
  # # create folds.json file
  # with open('../data/folds.json', 'w', encoding='utf8') as fp:
  #   json.dump(variants_folds, fp, indent=2, ensure_ascii=False)